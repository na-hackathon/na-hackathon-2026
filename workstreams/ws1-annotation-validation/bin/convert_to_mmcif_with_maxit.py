#!/usr/bin/env python3
"""
convert_to_mmcif_with_maxit.py

Convert a structure file to mmCIF using a pre-built MAXIT Docker image.
The input format is determined from the file extension.

Supported input formats:
    .pdb            -> converted to mmCIF
    .cif / .mmcif   -> run through maxit to make sure it's updated mmCIF

Usage:
    python pdb_to_mmcif.py <input_file> --image <image> [-o OUTPUT] [--force]

If --output is omitted, the output file is written in the same directory as
the input file with a .cif extension.

If the input is already mmCIF, the script reports this and exits successfully.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

PDB_EXTS = {".pdb"}
MMCIF_EXTS = {".cif", ".mmcif"}


def log(msg: str) -> None:
    print(f"[pdb_to_mmcif] {msg}", file=sys.stderr)


def die(msg: str) -> "SystemExit":
    return SystemExit(f"error: {msg}")


def convert_pdb_to_mmcif(input_path: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log(f"converting {input_path} -> {output_path}")
    
    # no docker call
    proc = subprocess.run(
        ['maxit', '-input', input_path, '-output', output_path, '-o', '1' ],
    )

    if proc.returncode != 0:
        output_path.unlink(missing_ok=True)  # don't leave a half-written file behind
        sys.stderr.write(proc.stderr.decode("utf-8", errors="replace"))
        raise die(f"maxit conversion failed (exit {proc.returncode})")

def convert_cif_to_mmcif(input_path: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log(f"converting {input_path} -> {output_path}")
    
    # no docker call
    proc = subprocess.run(
        ['maxit', '-input', input_path, '-output', output_path, '-o', '8' ],
    )

    if proc.returncode != 0:
        output_path.unlink(missing_ok=True)  # don't leave a half-written file behind
        sys.stderr.write(proc.stderr.decode("utf-8", errors="replace"))
        raise die(f"maxit conversion failed (exit {proc.returncode})")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert a PDB (or cif) file to (updated)mmCIF via a pre-built maxit from conda."
    )
    p.add_argument("input", type=Path, help="Path to input structure file.")
    p.add_argument("-o", "--output", type=Path,
                   help="Output path (default: <input>.cif next to the input).")
    p.add_argument("--force", action="store_true",
                   help="Overwrite the output file if it already exists.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    input_path: Path = args.input.resolve()
    if not input_path.is_file():
        raise die(f"input file not found: {input_path}")

    output_path: Path = (
        args.output.resolve() if args.output is not None
        else input_path.with_suffix(".cif")
    )

    ext = input_path.suffix.lower()
    if ext in MMCIF_EXTS:
        log(f"input is already CIF ({ext}) - let's update it")
        convert_cif_to_mmcif(input_path, output_path)
        return 0

    if ext not in PDB_EXTS:
        raise die(f"unrecognized extension {ext!r} for {input_path} "
                  f"(expected one of {sorted(PDB_EXTS | MMCIF_EXTS)})")

    convert_pdb_to_mmcif(input_path, output_path)
    log(f"done: {output_path}")
    print(output_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
