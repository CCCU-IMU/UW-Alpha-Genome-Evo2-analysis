#!/usr/bin/env python3
"""Scan simple PWMs around REF/ALT sequences and report motif gain/loss.

This lightweight scanner is intended for transparent examples. For manuscript
work, use a curated motif database such as JASPAR and document the motif release.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict

from uw_variant_prioritization.utils import read_tsv, to_float, write_tsv


BASES = ["A", "C", "G", "T"]


def read_fasta(path: str) -> dict[str, str]:
    records: dict[str, str] = {}
    name = ""
    seq: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records[name] = "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            records[name] = "".join(seq).upper()
    return records


def read_pwms(path: str) -> dict[str, list[dict[str, float]]]:
    rows = read_tsv(path)
    motifs: dict[str, list[dict[str, float]]] = defaultdict(list)
    for row in rows:
        motifs[row["motif_id"]].append({base: max(to_float(row.get(base, ""), 0.0), 1e-6) for base in BASES})
    return dict(motifs)


def score_window(seq: str, pwm: list[dict[str, float]]) -> float:
    score = 0.0
    for base, probs in zip(seq, pwm):
        score += math.log2(probs.get(base, 1e-6) / 0.25)
    return score


def best_score(seq: str, pwm: list[dict[str, float]]) -> float:
    width = len(pwm)
    if len(seq) < width:
        return math.nan
    return max(score_window(seq[i : i + width], pwm) for i in range(len(seq) - width + 1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", default="examples/ref_alt_sequences.example.fa", help="FASTA with IDs like variant_id|REF and variant_id|ALT.")
    parser.add_argument("--motifs", default="examples/motif_pwm.example.tsv", help="Simple PWM TSV with motif_id, position, A, C, G, T.")
    parser.add_argument("--output", default="results/targeted_motif_gain_loss.tsv", help="Output motif gain/loss table.")
    parser.add_argument("--dry-run", action="store_true", help="Read inputs and report planned output without writing.")
    args = parser.parse_args()

    sequences = read_fasta(args.fasta)
    motifs = read_pwms(args.motifs)
    variants = sorted({name.rsplit("|", 1)[0] for name in sequences})
    out = []
    for variant in variants:
        ref = sequences.get(f"{variant}|REF")
        alt = sequences.get(f"{variant}|ALT")
        if not ref or not alt:
            continue
        for motif_id, pwm in motifs.items():
            ref_score = best_score(ref, pwm)
            alt_score = best_score(alt, pwm)
            delta = alt_score - ref_score if math.isfinite(ref_score) and math.isfinite(alt_score) else math.nan
            out.append(
                {
                    "variant_id": variant,
                    "motif_id": motif_id,
                    "ref_score": f"{ref_score:.6g}",
                    "alt_score": f"{alt_score:.6g}",
                    "delta_ALT_minus_REF": f"{delta:.6g}",
                    "effect": "gain" if math.isfinite(delta) and delta > 0 else "loss" if math.isfinite(delta) and delta < 0 else "near-zero",
                }
            )
    if args.dry_run:
        print(f"Would scan {len(variants)} variants across {len(motifs)} motifs and write {args.output}")
        return
    fields = ["variant_id", "motif_id", "ref_score", "alt_score", "delta_ALT_minus_REF", "effect"]
    write_tsv(args.output, out, fields)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
