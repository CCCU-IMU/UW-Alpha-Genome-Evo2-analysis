#!/usr/bin/env python3
"""Infer UW-enriched allele direction from allele-frequency summaries.

The input is a small tabular allele-frequency summary, not a private VCF.
If using a VCF, first derive group AC/AN/AF values with your own population map.
"""

from __future__ import annotations

import argparse
import math

from uw_variant_prioritization.utils import read_tsv, to_float, write_tsv


def classify(row: dict[str, str], margin: float) -> tuple[str, str]:
    uw = to_float(row.get("UW_AF", ""), math.nan)
    controls = [
        to_float(row.get(field, ""), math.nan)
        for field in ["Angus_AF", "Charolais_AF", "MoOD_AF", "MoSN_AF", "Mongolia_AF", "control_AF"]
    ]
    controls = [x for x in controls if math.isfinite(x)]
    if not math.isfinite(uw) or not controls:
        return "low_confidence", "low_confidence"
    control = sum(controls) / len(controls)
    if uw - control >= margin:
        return "ALT", "UW_enriched_ALT"
    if control - uw >= margin:
        return "REF", "UW_enriched_REF"
    return "intermediate", "intermediate"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="examples/allele_frequency.example.tsv", help="Allele-frequency summary TSV.")
    parser.add_argument("--output", default="results/allele_direction.tsv", help="Output allele-direction TSV.")
    parser.add_argument("--margin", type=float, default=0.05, help="Minimum UW-control AF difference for direction calls.")
    parser.add_argument("--dry-run", action="store_true", help="Read inputs and report planned output without writing.")
    args = parser.parse_args()

    rows = read_tsv(args.input)
    out = []
    for row in rows:
        allele, direction = classify(row, args.margin)
        record = dict(row)
        record["UW_enriched_allele"] = allele
        record["direction_class"] = direction
        out.append(record)
    if args.dry_run:
        print(f"Would classify {len(out)} variants and write {args.output}")
        return
    fields = list(out[0].keys()) if out else []
    write_tsv(args.output, out, fields)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
