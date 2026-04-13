#!/usr/bin/env python3
"""Suggest next (amplitude, temperature) from run history for M0/M2 targeting.

Output format (stdout):
  <next_amp> <next_temp_eV> <method>
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np


@dataclass
class HistoryRow:
    amp: float
    temp_eV: float
    m0: float
    m2: float


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Suggest next source parameters")
    p.add_argument("--history", required=True, help="CSV history file")
    p.add_argument("--m0-target", type=float, required=True)
    p.add_argument("--m2-target", type=float, required=True)
    p.add_argument("--current-amp", type=float, required=True)
    p.add_argument("--current-temp", type=float, required=True)
    p.add_argument("--min-amp", type=float, default=1e-30)
    p.add_argument("--min-temp", type=float, default=1e-6)
    return p.parse_args()


def load_history(path: str) -> List[HistoryRow]:
    rows: List[HistoryRow] = []
    try:
        with open(path, "r", encoding="utf-8") as f:
            rdr = csv.DictReader(f)
            for r in rdr:
                try:
                    amp = float(r["amp"])
                    temp = float(r["temp_eV"])
                    m0 = float(r["m0"])
                    m2 = float(r["m2"])
                except (KeyError, ValueError):
                    continue
                if amp > 0 and temp > 0 and m0 > 0 and m2 > 0:
                    rows.append(HistoryRow(amp, temp, m0, m2))
    except FileNotFoundError:
        return []
    return rows


def heuristic_update(
    current_amp: float,
    current_temp: float,
    m0_target: float,
    m2_target: float,
    last_m0: float,
    last_m2: float,
) -> Tuple[float, float, str]:
    amp = current_amp
    temp = current_temp
    if last_m0 > 0:
        amp = current_amp * (m0_target / last_m0)
    if last_m2 > 0:
        temp = current_temp * (m2_target / last_m2)
    return amp, temp, "heuristic"


def fit_log_surrogate(rows: List[HistoryRow]):
    """Fit linear surrogate in log-space.

    log(M0) = a0 + a1*log(amp) + a2*log(temp)
    log(M2) = b0 + b1*log(amp) + b2*log(temp)
    """
    X = []
    y0 = []
    y2 = []
    for r in rows:
        X.append([1.0, math.log(r.amp), math.log(r.temp_eV)])
        y0.append(math.log(r.m0))
        y2.append(math.log(r.m2))

    Xn = np.array(X)
    y0n = np.array(y0)
    y2n = np.array(y2)

    a, *_ = np.linalg.lstsq(Xn, y0n, rcond=None)
    b, *_ = np.linalg.lstsq(Xn, y2n, rcond=None)
    return a, b


def suggest_with_scipy(
    rows: List[HistoryRow],
    current_amp: float,
    current_temp: float,
    m0_target: float,
    m2_target: float,
) -> Tuple[float, float, str]:
    try:
        from scipy.optimize import least_squares
    except Exception:
        return current_amp, current_temp, "no-scipy"

    a, b = fit_log_surrogate(rows)

    log_t0 = math.log(m0_target)
    log_t2 = math.log(m2_target)

    def residuals(logx: np.ndarray) -> np.ndarray:
        la, lt = float(logx[0]), float(logx[1])
        p0 = a[0] + a[1] * la + a[2] * lt
        p2 = b[0] + b[1] * la + b[2] * lt
        return np.array([p0 - log_t0, p2 - log_t2])

    amps = np.array([r.amp for r in rows])
    temps = np.array([r.temp_eV for r in rows])

    lb_amp = max(float(np.min(amps)) / 3.0, 1e-30)
    ub_amp = max(float(np.max(amps)) * 3.0, lb_amp * 1.01)
    lb_tmp = max(float(np.min(temps)) / 3.0, 1e-12)
    ub_tmp = max(float(np.max(temps)) * 3.0, lb_tmp * 1.01)

    x0 = np.array([math.log(max(current_amp, lb_amp)), math.log(max(current_temp, lb_tmp))])
    lower = np.array([math.log(lb_amp), math.log(lb_tmp)])
    upper = np.array([math.log(ub_amp), math.log(ub_tmp)])

    res = least_squares(residuals, x0=x0, bounds=(lower, upper), method="trf")
    amp = math.exp(float(res.x[0]))
    temp = math.exp(float(res.x[1]))
    return amp, temp, f"scipy-trf(cost={res.cost:.3e})"


def main() -> int:
    args = parse_args()

    rows = load_history(args.history)
    if not rows:
        print(f"{args.current_amp:.12g} {args.current_temp:.12g} keep-current")
        return 0

    last = rows[-1]

    # Use a simple heuristic for the first two samples, then switch to model-based search.
    if len(rows) < 3:
        amp, temp, method = heuristic_update(
            args.current_amp,
            args.current_temp,
            args.m0_target,
            args.m2_target,
            last.m0,
            last.m2,
        )
    else:
        amp, temp, method = suggest_with_scipy(
            rows,
            args.current_amp,
            args.current_temp,
            args.m0_target,
            args.m2_target,
        )
        if method == "no-scipy":
            amp, temp, method = heuristic_update(
                args.current_amp,
                args.current_temp,
                args.m0_target,
                args.m2_target,
                last.m0,
                last.m2,
            )
            method = "heuristic(no-scipy)"

    amp = max(float(amp), float(args.min_amp))
    temp = max(float(temp), float(args.min_temp))

    if not math.isfinite(amp) or not math.isfinite(temp):
        print(f"{args.current_amp:.12g} {args.current_temp:.12g} fallback-nonfinite")
        return 0

    print(f"{amp:.12g} {temp:.12g} {method}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
