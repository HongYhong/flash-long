"""Visualization CLI entrypoint (`isoquant-next-visualize`)."""

from __future__ import annotations

import sys
from typing import List, Optional


def main(cmd_args: Optional[List[str]] = None) -> int:
    from isoquant_next.visualization.cli import run_visualize

    argv = list(sys.argv[1:] if cmd_args is None else cmd_args)
    return int(run_visualize(argv))


if __name__ == "__main__":
    raise SystemExit(main())
