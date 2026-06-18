#!/usr/bin/env python3
"""Compute LD r^2 between representative SNPs and top delta-AF SNPs.

Inputs are lightweight dosage tables for demonstration. In a real analysis,
generate dosage values from your own VCF after applying your population filters.
"""

from __future__ import annotations

import argparse
import math

from uw_variant_prioritization.utils import read_tsv, to_float, write_tsv


def r2(xs: list[float], ys: list[float]) -> float:
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 3:
        return math.nan
    xvals = [p[0] for p in pairs]
    yvals = [p[1] for p in pairs]
    xmean = sum(xvals) / len(xvals)
    ymean = sum(yvals) / len(yvals)
    cov = sum((x - xmean) * (y - ymean) for x, y in pairs)
    vx = sum((x - xmean) ** 2 for x in xvals)
    vy = sum((y - ymean) ** 2 for y in yvals)
    if vx == 0 or vy == 0:
        return math.nan
    return (cov * cov) / (vx * vy)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", default="examples/ld_pairs.example.tsv", help="Representative/top-variant pairs TSV.")
    parser.add_argument("--dosage", default="examples/dosage.example.tsv", help="Wide dosage table, variant_id plus sample columns.")
    parser.add_argument("--output", default="results/representative_vs_top_deltaAF_LD.tsv", help="Output LD table.")
    parser.add_argument("--threshold", type=float, default=0.30, help="r^2 threshold for LD-haplotype support.")
    parser.add_argument("--dry-run", action="store_true", help="Read inputs and report planned output without writing.")
    args = parser.parse_args()

    pairs = read_tsv(args.pairs)
    dosage_rows = read_tsv(args.dosage)
    dosage = {row["variant_id"]: row for row in dosage_rows}
    sample_cols = [c for c in dosage_rows[0] if c != "variant_id"] if dosage_rows else []
    out = []
    for pair in pairs:
        rep = dosage.get(pair["representative_variant"])
        top = dosage.get(pair["top_deltaAF_variant"])
        value = math.nan
        if rep and top:
            value = r2([to_float(rep.get(s, ""), math.nan) for s in sample_cols], [to_float(top.get(s, ""), math.nan) for s in sample_cols])
        record = dict(pair)
        record["ld_r2"] = "" if not math.isfinite(value) else f"{value:.6g}"
        record["support_class"] = "LD-haplotype" if math.isfinite(value) and value >= args.threshold else "same-window_or_low_LD"
        out.append(record)
    if args.dry_run:
        print(f"Would compute LD for {len(out)} pairs and write {args.output}")
        return
    fields = list(out[0].keys()) if out else []
    write_tsv(args.output, out, fields)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
