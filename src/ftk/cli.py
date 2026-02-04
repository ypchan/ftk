from __future__ import annotations

import argparse
import sys

from rich_argparse import RichHelpFormatter

from ftk.commands.stat import add_stat_parser


def _build_parser() -> argparse.ArgumentParser:
    # RichHelpFormatter will render argparse help via Rich.
    # add_argument(..., default=...) will be shown automatically by ArgumentDefaultsHelpFormatter behavior.
    # RichHelpFormatter already supports showing defaults; keep defaults explicit on args.
    parser = argparse.ArgumentParser(
        prog="ftk",
        description="ftk: extensible CLI toolkit.",
        formatter_class=RichHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest="command", required=True)
    add_stat_parser(subparsers)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "stat":
        from ftk.commands.stat import run_stat
        run_stat(args)
        return

    parser.print_help()


if __name__ == "__main__":
    main(sys.argv[1:])
