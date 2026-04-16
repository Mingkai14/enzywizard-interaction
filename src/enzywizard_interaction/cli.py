from __future__ import annotations

import argparse

from .commands.interaction import add_interaction_parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="enzywizard-interaction",
        description="EnzyWizard-Interaction: Calculate protein/protein-substrate interactions and generate a detailed JSON report."
    )
    add_interaction_parser(parser)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)