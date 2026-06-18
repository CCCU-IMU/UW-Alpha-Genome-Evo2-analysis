"""Shared utility functions for the public UW prioritization workflow."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable


def readable_path(path: str | Path) -> Path:
    p = Path(path)
    if not p.exists():
        raise argparse.ArgumentTypeError(f"file not found: {p}")
    return p


def ensure_parent(path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: str | Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    p = ensure_parent(path)
    with p.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def to_float(value: object, default: float = math.nan) -> float:
    if value is None:
        return default
    text = str(value).strip()
    if text == "" or text.upper() in {"NA", "NAN", "NONE"}:
        return default
    try:
        return float(text)
    except ValueError:
        return default


def minmax(values: list[float]) -> list[float]:
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return [math.nan for _ in values]
    lo = min(finite)
    hi = max(finite)
    if hi == lo:
        return [1.0 if math.isfinite(v) else math.nan for v in values]
    return [(v - lo) / (hi - lo) if math.isfinite(v) else math.nan for v in values]


def log1p_minmax(values: list[float]) -> list[float]:
    transformed = [math.log1p(max(v, 0.0)) if math.isfinite(v) else math.nan for v in values]
    return minmax(transformed)


def weighted_mean(values: list[float], weights: list[float]) -> float:
    used = [(v, w) for v, w in zip(values, weights) if math.isfinite(v) and w > 0]
    if not used:
        return math.nan
    return sum(v * w for v, w in used) / sum(w for _v, w in used)


def effect_for_uw_allele(alt_minus_ref: float, uw_allele: str) -> float:
    """Convert an ALT-REF model direction to the UW-biased allele direction."""
    allele = (uw_allele or "").upper()
    if not math.isfinite(alt_minus_ref):
        return math.nan
    if allele == "ALT":
        return alt_minus_ref
    if allele == "REF":
        return -alt_minus_ref
    return math.nan


def effect_label(value: float, threshold: float = 1e-9) -> str:
    if not math.isfinite(value):
        return "unresolved"
    if value > threshold:
        return "increase"
    if value < -threshold:
        return "decrease"
    return "near-zero"
