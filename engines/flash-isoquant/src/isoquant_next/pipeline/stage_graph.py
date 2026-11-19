"""Declarative Stage Graph DAG with atomic JSON state tracking for deterministic --resume."""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import logging
import os
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from isoquant_next.core.error_codes import IsoQuantError, IsoQuantExitCode

logger = logging.getLogger("IsoQuant")


@dataclass(frozen=True, slots=True)
class StageNode:
    """Single pipeline stage node in the topological DAG."""

    name: str
    run: Callable[[], None]
    deps: Tuple[str, ...] = ()
    enabled: bool = True
    on_skip: Optional[Callable[[], None]] = None


class StageStateTracker:
    """Persists completed pipeline stages atomically via os.replace for safe --resume."""

    def __init__(self, state_path: str | Path, fingerprint: str) -> None:
        self.path = Path(state_path)
        self.fingerprint = fingerprint
        self.completed: Dict[str, dict] = {}

    def load(self) -> None:
        if not self.path.exists():
            self.completed = {}
            return
        try:
            doc = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise IsoQuantError(
                IsoQuantExitCode.FILE_FORMAT_ERROR,
                f"Corrupted stage state file {self.path}: {exc}",
            ) from exc
        if doc.get("fingerprint") != self.fingerprint:
            logger.info("Pipeline configuration fingerprint changed; resetting stage completion state")
            self.completed = {}
            return
        self.completed = dict(doc.get("completed_stages", {}))

    def is_done(self, stage_name: str) -> bool:
        return stage_name in self.completed

    def mark_done(self, stage_name: str, metadata: Optional[dict] = None) -> None:
        self.completed[stage_name] = metadata or {"status": "done"}
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(self.path.suffix + f".tmp.{os.getpid()}")
        payload = {
            "fingerprint": self.fingerprint,
            "completed_stages": self.completed,
        }
        tmp.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        os.replace(tmp, self.path)

    def clear(self) -> None:
        self.completed.clear()
        if self.path.exists():
            self.path.unlink()


class StageGraphExecutor:
    """Executes enabled StageNodes in topological order with resume skipping and cycle validation."""

    def __init__(self, stages: List[StageNode], tracker: Optional[StageStateTracker] = None) -> None:
        self.stages = stages
        self.tracker = tracker
        self._order = self._toposort(stages)

    @staticmethod
    def _toposort(stages: List[StageNode]) -> List[StageNode]:
        by_name = {s.name: s for s in stages}
        indeg: Dict[str, int] = {s.name: 0 for s in stages}
        children: Dict[str, List[str]] = {s.name: [] for s in stages}
        for s in stages:
            for dep in s.deps:
                if dep not in by_name:
                    raise ValueError(f"Unknown stage dependency '{dep}' for stage '{s.name}'")
                indeg[s.name] += 1
                children[dep].append(s.name)
        queue = [s.name for s in stages if indeg[s.name] == 0]
        ordered: List[StageNode] = []
        while queue:
            cur = queue.pop(0)
            ordered.append(by_name[cur])
            for nxt in children[cur]:
                indeg[nxt] -= 1
                if indeg[nxt] == 0:
                    queue.append(nxt)
        if len(ordered) != len(stages):
            raise ValueError("Cycle detected in pipeline StageGraph")
        return ordered

    def execute(self, *, resume: bool = False) -> List[str]:
        executed: List[str] = []
        if self.tracker is not None:
            if resume:
                self.tracker.load()
            else:
                self.tracker.clear()

        for node in self._order:
            if not node.enabled:
                continue
            if resume and self.tracker is not None and self.tracker.is_done(node.name):
                logger.info("Resuming: skipping completed stage '%s'", node.name)
                if node.on_skip is not None:
                    node.on_skip()
                continue
            node.run()
            executed.append(node.name)
            if self.tracker is not None:
                self.tracker.mark_done(node.name)
        return executed
