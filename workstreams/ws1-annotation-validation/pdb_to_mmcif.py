#!/usr/bin/env python3
"""
pdb_to_mmcif.py - given a path to a structure file, convert it to mmCIF using a
pre-built maxit Docker image. Format is determined by file extension:

    .pdb / .ent            -> PDB,   converted to mmCIF
    .cif / .mmcif          -> mmCIF, left as-is (no-op)

Usage:
    ./pdb_to_mmcif.py <input_file> --image <image> [-o OUTPUT] [--force]

If --output is omitted, the result is written next to the input with a .cif
extension. If the input is already mmCIF, the script reports that and exits 0.
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


def convert_pdb_to_mmcif(input_path: Path, output_path: Path, image: str) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log(f"converting {input_path} -> {output_path}")

    with input_path.open("rb") as fin, output_path.open("wb") as fout:
        proc = subprocess.run(
            ["docker", "run", "--rm", "-i", image],
            stdin=fin, stdout=fout, stderr=subprocess.PIPE,
        )

    if proc.returncode != 0:
        output_path.unlink(missing_ok=True)  # don't leave a half-written file behind
        sys.stderr.write(proc.stderr.decode("utf-8", errors="replace"))
        raise die(f"maxit conversion failed (exit {proc.returncode})")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert a PDB file to mmCIF via a pre-built maxit Docker "
                    "image. Format is detected from the file extension.",
    )
    p.add_argument("input", type=Path, help="Path to input structure file.")
    p.add_argument("--image", default=os.environ.get("MAXIT_IMAGE"),
                   help="Docker image tag/ID to use (or set MAXIT_IMAGE env var).")
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

    ext = input_path.suffix.lower()
    if ext in MMCIF_EXTS:
        log(f"input is already mmCIF ({ext}) - nothing to do.")
        return 0
    if ext not in PDB_EXTS:
        raise die(f"unrecognized extension {ext!r} for {input_path} "
                  f"(expected one of {sorted(PDB_EXTS | MMCIF_EXTS)})")

    output_path: Path = (
        args.output.resolve() if args.output is not None
        else input_path.with_suffix(".cif")
    )

    if not args.image:
        raise die("no Docker image specified (pass --image or set MAXIT_IMAGE)")
    if output_path.exists() and not args.force:
        raise die(f"output already exists: {output_path} (use --force to overwrite)")
    if shutil.which("docker") is None:
        raise die("'docker' not found in PATH")

    convert_pdb_to_mmcif(input_path, output_path, args.image)
    log(f"done: {output_path}")
    print(output_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
