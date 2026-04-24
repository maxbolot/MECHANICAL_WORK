import os
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Optional

import numpy as np


@dataclass(frozen=True)
class TimeWeightResult:
    weights_days: np.ndarray
    missing_steps: int


def weighted_nanmean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=float).reshape(-1)
    weights = np.asarray(weights, dtype=float).reshape(-1)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(valid):
        return float("nan")
    return float(np.sum(values[valid] * weights[valid]) / np.sum(weights[valid]))


def normalize_time_axis_to_days(time_axis: np.ndarray, units: Optional[str]) -> np.ndarray:
    t = np.asarray(time_axis, dtype=float).reshape(-1)
    if units and "since" in units.lower():
        unit_token, ref = units.split("since", 1)
        unit_token = unit_token.strip().lower()
        ref_dt = _parse_reference_time(ref.strip())
        ref_days = ref_dt.timestamp() / 86400.0
        if unit_token in {"days", "day"}:
            return ref_days + t
        if unit_token in {"hours", "hour", "hrs", "hr"}:
            return ref_days + t / 24.0
        if unit_token in {"minutes", "minute", "mins", "min"}:
            return ref_days + t / (24.0 * 60.0)
        if unit_token in {"seconds", "second", "secs", "sec"}:
            return ref_days + t / (24.0 * 60.0 * 60.0)

    if np.all((t > 1.0e6) & (t < 1.0e10)):
        day_part = np.floor(t).astype(int)
        frac_part = t - day_part
        out = np.empty_like(t, dtype=float)
        for i, yyyymmdd in enumerate(day_part):
            dt = datetime.strptime(str(yyyymmdd), "%Y%m%d").replace(tzinfo=timezone.utc)
            out[i] = dt.timestamp() / 86400.0 + frac_part[i]
        return out

    raise ValueError("Unable to convert time axis to day numbers")


def compute_time_weights_from_days(time_days: np.ndarray, tag: str = "simulation") -> TimeWeightResult:
    td = np.asarray(time_days, dtype=float).reshape(-1)
    ntime = td.size
    if ntime == 0:
        return TimeWeightResult(weights_days=np.array([], dtype=float), missing_steps=0)
    if not np.all(np.isfinite(td)):
        raise ValueError(f"Invalid non-finite time values in {tag} time axis")

    td = np.sort(td)
    if ntime > 1:
        dt = np.diff(td)
        if np.any(dt <= 0):
            raise ValueError(f"Time axis is not strictly increasing ({tag})")
        if np.median(dt) > 30:
            raise ValueError(f"Unrealistic median timestep ({np.median(dt)} days) in {tag}")

    tol_hours = _get_missing_timestep_tolerance_hours()
    tol_days = tol_hours / 24.0
    missing_steps = _count_missing_daily_schedule(td, tol_days)
    return TimeWeightResult(weights_days=np.ones(ntime, dtype=float), missing_steps=missing_steps)


def compute_time_weights_control(time_axis: np.ndarray, units: Optional[str]) -> TimeWeightResult:
    return compute_time_weights_from_days(normalize_time_axis_to_days(time_axis, units), tag="control")


def compute_time_weights_plus4k(time_axis: np.ndarray, units: Optional[str]) -> TimeWeightResult:
    return compute_time_weights_from_days(normalize_time_axis_to_days(time_axis, units), tag="plus4k")


def _count_missing_daily_schedule(time_days: np.ndarray, tol_days: float) -> int:
    if time_days.size <= 1:
        return 0
    expected = np.arange(time_days[0], time_days[-1] + 1.0e-12, 1.0)
    missing_steps = 0
    iobs = 0
    nobs = time_days.size
    for exp in expected:
        while iobs < nobs and time_days[iobs] < exp - tol_days:
            iobs += 1
        if iobs >= nobs or abs(time_days[iobs] - exp) > tol_days:
            missing_steps += 1
    return int(missing_steps)


def _get_missing_timestep_tolerance_hours() -> float:
    value = os.getenv("MISSING_TIMESTEP_TOL_HOURS", "18")
    tol_hours = float(value)
    if tol_hours <= 0:
        raise ValueError("MISSING_TIMESTEP_TOL_HOURS must be positive")
    return tol_hours


def _parse_reference_time(value: str) -> datetime:
    text = value.strip().replace("UTC", "").replace("Z", "")
    fmts = [
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d %H:%M",
        "%Y-%m-%d",
        "%Y/%m/%d %H:%M:%S",
        "%Y/%m/%d %H:%M",
        "%Y/%m/%d",
    ]
    for fmt in fmts:
        try:
            return datetime.strptime(text, fmt).replace(tzinfo=timezone.utc)
        except ValueError:
            continue
    raise ValueError(f"Unsupported reference time format: {value}")
