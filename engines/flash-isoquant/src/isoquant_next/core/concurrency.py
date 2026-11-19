"""Safe multiprocessing runner abstracting sequential / forkserver / spawn pools."""

from __future__ import annotations

import logging
import multiprocessing as mp
import sys
import tempfile
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from typing import Callable, Iterable, Iterator, List, Optional, Sequence, TypeVar

from isoquant_next.core.logging_setup import configure_logging, pin_single_thread_math_env

T = TypeVar("T")
R = TypeVar("R")


def _worker_bootstrap(debug: bool, log_file: Optional[str], custom_init: Optional[Callable[[], None]]) -> None:
    pin_single_thread_math_env(force=True)
    configure_logging(debug=debug, log_file=log_file)
    if custom_init is not None:
        custom_init()


def select_mp_context(preferred: str = "forkserver") -> str:
    """Choose a forkserver/spawn context safe with HTSlib and BLAS native libraries."""
    available = mp.get_all_start_methods()
    if preferred == "forkserver" and "forkserver" in available and sys.platform.startswith("linux"):
        # Linux AF_UNIX listener paths inside forkserver break when TMPDIR exceeds ~65 chars.
        if len(tempfile.gettempdir()) <= 65:
            return "forkserver"
        if "spawn" in available:
            return "spawn"
    if preferred in available:
        return preferred
    return "spawn" if "spawn" in available else available[0]


class WorkerPoolExecutor:
    """Unified task executor replacing bare ProcessPoolExecutor(mp_context=fork) scatter."""

    def __init__(
        self,
        max_workers: int = 1,
        mp_method: str = "forkserver",
        debug: bool = False,
        log_file: Optional[str] = None,
        initializer: Optional[Callable[[], None]] = None,
    ) -> None:
        self.max_workers = max(1, int(max_workers))
        self.mp_method = mp_method
        self.debug = debug
        self.log_file = log_file
        self.initializer = initializer
        self._pool: Optional[ProcessPoolExecutor] = None

    def __enter__(self) -> "WorkerPoolExecutor":
        pin_single_thread_math_env(force=True)
        if self.max_workers > 1:
            ctx_name = select_mp_context(self.mp_method)
            ctx = mp.get_context(ctx_name)
            self._pool = ProcessPoolExecutor(
                max_workers=self.max_workers,
                mp_context=ctx,
                initializer=_worker_bootstrap,
                initargs=(self.debug, self.log_file, self.initializer),
            )
        elif self.initializer is not None:
            self.initializer()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        if self._pool is not None:
            self._pool.shutdown(wait=True)
            self._pool = None

    def submit(self, fn: Callable[..., R], *args, **kwargs) -> Future[R]:
        if self._pool is not None:
            return self._pool.submit(fn, *args, **kwargs)
        fut: Future[R] = Future()
        try:
            res = fn(*args, **kwargs)
        except BaseException as exc:
            fut.set_exception(exc)
        else:
            fut.set_result(res)
        return fut

    def map_ordered(self, fn: Callable[[T], R], items: Sequence[T]) -> List[R]:
        if not items:
            return []
        if self.max_workers <= 1 or len(items) == 1:
            pin_single_thread_math_env(force=True)
            if self.initializer is not None:
                self.initializer()
            return [fn(item) for item in items]

        ctx_name = select_mp_context(self.mp_method)
        ctx = mp.get_context(ctx_name)
        pin_single_thread_math_env(force=True)
        with ProcessPoolExecutor(
            max_workers=min(self.max_workers, len(items)),
            mp_context=ctx,
            initializer=_worker_bootstrap,
            initargs=(self.debug, self.log_file, self.initializer),
        ) as pool:
            return list(pool.map(fn, items))

    def map_unordered(self, fn: Callable[[T], R], items: Iterable[T]) -> Iterator[R]:
        items_list = list(items)
        if not items_list:
            return
        if self.max_workers <= 1 or len(items_list) == 1:
            pin_single_thread_math_env(force=True)
            if self.initializer is not None:
                self.initializer()
            for item in items_list:
                yield fn(item)
            return

        ctx_name = select_mp_context(self.mp_method)
        ctx = mp.get_context(ctx_name)
        pin_single_thread_math_env(force=True)
        with ProcessPoolExecutor(
            max_workers=min(self.max_workers, len(items_list)),
            mp_context=ctx,
            initializer=_worker_bootstrap,
            initargs=(self.debug, self.log_file, self.initializer),
        ) as pool:
            futures = [pool.submit(fn, item) for item in items_list]
            try:
                for fut in as_completed(futures):
                    yield fut.result()
            except Exception:
                for f in futures:
                    f.cancel()
                logging.getLogger("IsoQuant").exception("Worker pool execution failed")
                raise
