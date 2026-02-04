from __future__ import annotations

import argparse
import sys

from rich_argparse import RichHelpFormatter

from ftk.commands.stat import add_stat_parser
from ftk.commands.fishing import add_fishing_parser
from ftk.commands.msaqc import add_msaqc_parser


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ftk",
        description="ftk: extensible CLI toolkit.",
        formatter_class=RichHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    add_stat_parser(subparsers)
    add_fishing_parser(subparsers)
    add_msaqc_parser(subparsers)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "stat":
        from ftk.commands.stat import run_stat
        run_stat(args)
        return

    if args.command == "fishing":
        from ftk.commands.fishing import run_fishing
        run_fishing(args)
        return
    
    if args.command == "msaqc":
        from ftk.commands.msaqc import run_msaqc
        run_msaqc(args)
        return


    parser.print_help()


if __name__ == "__main__":
    main(sys.argv[1:])
