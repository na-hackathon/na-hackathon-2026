#!/usr/bin/env python3
"""
dnatco_fetch.py

Download DNATCO "extended" mmCIF files for a set of PDB ids.

URL layout (verified against dnatco.datmos.org):
    https://dnatco.datmos.org/db/coordinates/<pad[:4]>/<pdbid[1:3]>/<pad>_dnatco_extended.cif.gz
where  pad = pdbid.lower().zfill(8)   (e.g. "1ehz" -> "00001ehz", shard "0000"/"eh").

These files are full coordinate mmCIFs that already embed NDB base-pair
annotations (used downstream by dnatco_to_graph.py).

Get ids either from a file (one per line), the CLI, or directly from RCSB:
    python dnatco_fetch.py --rcsb-rna 1000 -o cifs/
    python dnatco_fetch.py --ids-file ids.txt -o cifs/
    python dnatco_fetch.py 1ehz 1fmn 2q1r -o cifs/
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

BASE_URL = "https://dnatco.datmos.org/db/coordinates"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"


def dnatco_url(pdbid: str) -> str:
    pid = pdbid.lower()
    pad = pid.zfill(8)
    return f"{BASE_URL}/{pad[:4]}/{pid[1:3]}/{pad}_dnatco_extended.cif.gz"


def rcsb_rna_ids(limit: int) -> list[str]:
    """Fetch PDB ids of entries containing at least one RNA polymer entity."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.polymer_entity_count_RNA",
                "operator": "greater",
                "value": 0,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit}},
    }
    url = RCSB_SEARCH + "?json=" + urllib.parse.quote(json.dumps(query))
    with urllib.request.urlopen(url, timeout=60) as r:
        data = json.load(r)
    return [hit["identifier"].lower() for hit in data.get("result_set", [])]


def fetch_one(pdbid: str, out_dir: Path, force: bool) -> tuple[str, str]:
    out = out_dir / f"{pdbid.lower()}.cif"
    if out.exists() and not force:
        return pdbid, "skip"
    try:
        with urllib.request.urlopen(dnatco_url(pdbid), timeout=120) as r:
            raw = r.read()
        out.write_bytes(gzip.decompress(raw))
        return pdbid, "ok"
    except urllib.error.HTTPError as e:
        return pdbid, f"http_{e.code}"
    except Exception as e:  # noqa: BLE001 - report and continue
        return pdbid, f"err_{type(e).__name__}"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("ids", nargs="*", help="PDB ids to fetch.")
    p.add_argument("--ids-file", type=Path, help="File with one PDB id per line.")
    p.add_argument("--rcsb-rna", type=int, metavar="N",
                   help="Fetch the first N RNA-containing PDB ids from RCSB.")
    p.add_argument("-o", "--out", type=Path, default=Path("cifs"),
                   help="Output directory for .cif files (default: cifs/).")
    p.add_argument("--workers", type=int, default=8, help="Concurrent downloads.")
    p.add_argument("--force", action="store_true", help="Re-download existing files.")
    return p.parse_args(argv)


def collect_ids(args: argparse.Namespace) -> list[str]:
    ids: list[str] = list(args.ids)
    if args.ids_file:
        ids += [ln.strip() for ln in args.ids_file.read_text().splitlines() if ln.strip()]
    if args.rcsb_rna:
        print(f"[fetch] querying RCSB for {args.rcsb_rna} RNA ids...", file=sys.stderr)
        ids += rcsb_rna_ids(args.rcsb_rna)
    # de-dup, keep order
    seen, uniq = set(), []
    for i in ids:
        i = i.lower()
        if i not in seen:
            seen.add(i)
            uniq.append(i)
    return uniq


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    ids = collect_ids(args)
    if not ids:
        raise SystemExit("error: no ids (pass ids, --ids-file, or --rcsb-rna).")
    args.out.mkdir(parents=True, exist_ok=True)

    counts: dict[str, int] = {}
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = [ex.submit(fetch_one, pid, args.out, args.force) for pid in ids]
        for n, fut in enumerate(cf.as_completed(futures), 1):
            pid, status = fut.result()
            counts[status.split("_")[0]] = counts.get(status.split("_")[0], 0) + 1
            if status.startswith(("err", "http")):
                print(f"[fetch] {pid}: {status}", file=sys.stderr)
            if n % 50 == 0:
                print(f"[fetch] {n}/{len(ids)} ...", file=sys.stderr)

    print(f"[fetch] done: {dict(counts)} -> {args.out}", file=sys.stderr)
    print(args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
