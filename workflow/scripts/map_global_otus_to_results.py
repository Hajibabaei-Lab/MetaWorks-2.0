#!/usr/bin/env python3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from map_global_to_results import main as _main


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    _main(argv + ["--column-name", "GlobalOTU"])


if __name__ == "__main__":
    main()
