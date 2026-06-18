#!/usr/bin/env python3
"""Create compact evidence-summary figures from an integrated evidence table."""

from __future__ import annotations

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np

from uw_variant_prioritization.utils import log1p_minmax, minmax, read_tsv, to_float


METRICS = [
    ("Selection", "priority_score", False),
    ("Evo2", "evo2_abs_delta_likelihood", False),
    ("DNASE", "dnase_center_abs_delta_sum", True),
    ("ATAC", "atac_center_abs_delta_sum", True),
    ("Expression", "expression_related_abs_score", True),
    ("Motif", "targeted_motif_abs_normalized_delta", True),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="examples/integrated_evidence_table.example.tsv", help="Integrated evidence table TSV.")
    parser.add_argument("--output-prefix", default="figures/example_evidence", help="Output prefix without extension.")
    parser.add_argument("--dry-run", action="store_true", help="Read inputs and report planned output without writing.")
    args = parser.parse_args()

    rows = read_tsv(args.input)
    if args.dry_run:
        print(f"Would plot {len(rows)} variants and write {args.output_prefix}.png/.pdf")
        return
    labels = [f"{r.get('candidate_gene', '')}\n{r.get('pos', '')}" for r in rows]
    matrix = []
    for _label, field, log_scale in METRICS:
        values = [to_float(r.get(field, ""), math.nan) for r in rows]
        matrix.append(log1p_minmax(values) if log_scale else minmax(values))
    arr = np.asarray(matrix, dtype=float)

    fig, ax = plt.subplots(figsize=(8.0, 3.6))
    im = ax.imshow(np.ma.masked_invalid(arr), aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_yticks(range(len(METRICS)))
    ax.set_yticklabels([m[0] for m in METRICS])
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax.tick_params(direction="in", length=0)
    cbar = fig.colorbar(im, ax=ax, pad=0.01)
    cbar.set_label("normalized evidence")
    ax.set_title("UW candidate variant model-based evidence matrix")
    fig.tight_layout()
    fig.savefig(f"{args.output_prefix}.png", dpi=300)
    fig.savefig(f"{args.output_prefix}.pdf")
    print(f"Wrote {args.output_prefix}.png and {args.output_prefix}.pdf")


if __name__ == "__main__":
    main()
