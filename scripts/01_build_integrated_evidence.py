#!/usr/bin/env python3
"""Build a normalized integrated evidence table for candidate variants.

This script expects precomputed model/output summaries rather than raw private
VCF, AlphaGenome API responses, or Evo2 model checkpoints. Long-tailed
regulatory metrics are log1p-transformed before min-max normalization.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from uw_variant_prioritization.utils import log1p_minmax, minmax, read_tsv, to_float, weighted_mean, write_tsv


DEFAULT_WEIGHTS = {
    "evo2": 0.40,
    "dnase": 0.20,
    "atac": 0.16,
    "expression": 0.06,
    "selection": 0.10,
    "motif": 0.08,
}


def parse_weights(text: str) -> dict[str, float]:
    weights = DEFAULT_WEIGHTS.copy()
    if not text:
        return weights
    for item in text.split(","):
        key, value = item.split("=", 1)
        weights[key.strip()] = float(value)
    return weights


def build_table(rows: list[dict[str, str]], weights: dict[str, float]) -> list[dict[str, object]]:
    metrics = {
        "selection": ("priority_score", minmax),
        "evo2": ("evo2_abs_delta_likelihood", minmax),
        "dnase": ("dnase_center_abs_delta_sum", log1p_minmax),
        "atac": ("atac_center_abs_delta_sum", log1p_minmax),
        "expression": ("expression_related_abs_score", log1p_minmax),
        "motif": ("targeted_motif_abs_normalized_delta", log1p_minmax),
    }
    normalized: dict[str, list[float]] = {}
    for key, (field, scaler) in metrics.items():
        normalized[key] = scaler([to_float(row.get(field, ""), math.nan) for row in rows])

    out: list[dict[str, object]] = []
    for i, row in enumerate(rows):
        vals = [normalized[key][i] for key in DEFAULT_WEIGHTS]
        wts = [weights[key] for key in DEFAULT_WEIGHTS]
        score = weighted_mean(vals, wts)
        record = dict(row)
        for key in DEFAULT_WEIGHTS:
            value = normalized[key][i]
            record[f"norm_{key}"] = "" if not math.isfinite(value) else f"{value:.6g}"
        record["integrated_score"] = "" if not math.isfinite(score) else f"{score:.6g}"
        out.append(record)
    out.sort(key=lambda r: to_float(r.get("integrated_score", ""), -1), reverse=True)
    for rank, row in enumerate(out, start=1):
        row["integrated_rank"] = rank
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="examples/integrated_evidence_table.example.tsv", help="Input evidence table TSV.")
    parser.add_argument("--output", default="results/integrated_evidence.ranked.tsv", help="Output ranked table TSV.")
    parser.add_argument("--weights", default="", help="Comma-separated weights, e.g. evo2=0.4,dnase=0.2.")
    parser.add_argument("--dry-run", action="store_true", help="Read inputs and report planned output without writing.")
    args = parser.parse_args()

    rows = read_tsv(args.input)
    ranked = build_table(rows, parse_weights(args.weights))
    if args.dry_run:
        print(f"Would rank {len(rows)} variants and write {args.output}")
        print("Required columns include priority_score, evo2_abs_delta_likelihood, DNASE/ATAC/expression/motif metrics.")
        return
    fields = list(ranked[0].keys()) if ranked else []
    write_tsv(args.output, ranked, fields)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
