from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from translated_from_matlab.lib.io import (
    ThresholdedCoordinates,
    read_required_var,
    read_standard_coordinates,
    read_thresholded_coordinates,
    read_thresholded_var,
)
from translated_from_matlab.lib.region_mask import RegionMask
from translated_from_matlab.lib.weights import weighted_nanmean


@dataclass(frozen=True)
class ThresholdedSpatialTimeseries:
    time: np.ndarray
    time_units: Optional[str]
    percentile: np.ndarray
    ntime: int
    nperc: int
    work_spatial_avg: np.ndarray
    lift_spatial_avg: np.ndarray


def compute_thresholded_region_spatial_timeseries(ncfile: str, mask: RegionMask) -> ThresholdedSpatialTimeseries:
    coords: ThresholdedCoordinates = read_thresholded_coordinates(ncfile)
    work4d = read_thresholded_var(ncfile, "work")
    lift4d = read_thresholded_var(ncfile, "lift")

    ntime = work4d.shape[2]
    nperc = work4d.shape[3]
    if nperc != coords.percentile.size:
        raise ValueError(
            f"Percentile dimension mismatch: data={nperc}, coordinate={coords.percentile.size}"
        )

    work_region = work4d[np.ix_(mask.lon_idx, mask.lat_idx, np.arange(ntime), np.arange(nperc))]
    lift_region = lift4d[np.ix_(mask.lon_idx, mask.lat_idx, np.arange(ntime), np.arange(nperc))]
    weight_4d = mask.lat_weights.reshape(1, -1, 1, 1)

    work_num = np.nansum(work_region * weight_4d, axis=(0, 1))
    work_den = np.nansum((~np.isnan(work_region)).astype(float) * weight_4d, axis=(0, 1))
    work_spatial_avg = work_num / work_den

    lift_num = np.nansum(lift_region * weight_4d, axis=(0, 1))
    lift_den = np.nansum((~np.isnan(lift_region)).astype(float) * weight_4d, axis=(0, 1))
    lift_spatial_avg = lift_num / lift_den

    return ThresholdedSpatialTimeseries(
        time=coords.time,
        time_units=coords.time_units,
        percentile=coords.percentile,
        ntime=ntime,
        nperc=nperc,
        work_spatial_avg=np.asarray(work_spatial_avg, dtype=float),
        lift_spatial_avg=np.asarray(lift_spatial_avg, dtype=float),
    )


def compute_unthresholded_region_time_means(cfg, mask: RegionMask) -> Tuple[float, float]:
    ncfile = cfg.scenario.standard_ncfile
    lat, lon, time, time_units = read_standard_coordinates(ncfile)

    work3d = read_required_var(ncfile, "work")
    lift3d = read_required_var(ncfile, "lift")

    work3d = _orient_to_lon_lat_time(work3d, nlon=lon.size, nlat=lat.size, ntime=time.size, varname="work")
    lift3d = _orient_to_lon_lat_time(lift3d, nlon=lon.size, nlat=lat.size, ntime=time.size, varname="lift")

    work_region = work3d[np.ix_(mask.lon_idx, mask.lat_idx, np.arange(work3d.shape[2]))]
    lift_region = lift3d[np.ix_(mask.lon_idx, mask.lat_idx, np.arange(lift3d.shape[2]))]
    weight_3d = mask.lat_weights.reshape(1, -1, 1)

    work_num = np.nansum(work_region * weight_3d, axis=(0, 1))
    work_den = np.nansum((~np.isnan(work_region)).astype(float) * weight_3d, axis=(0, 1))
    work_spatial_avg = work_num / work_den

    lift_num = np.nansum(lift_region * weight_3d, axis=(0, 1))
    lift_den = np.nansum((~np.isnan(lift_region)).astype(float) * weight_3d, axis=(0, 1))
    lift_spatial_avg = lift_num / lift_den

    time_weights = cfg.scenario.time_weight_fn(time, time_units).weights_days
    work_avg = weighted_nanmean(work_spatial_avg, time_weights)
    lift_avg = weighted_nanmean(lift_spatial_avg, time_weights)
    return work_avg, lift_avg


def _orient_to_lon_lat_time(data: np.ndarray, nlon: int, nlat: int, ntime: int, varname: str) -> np.ndarray:
    """Reorder a 3D field to (lon, lat, time) using coordinate lengths."""
    arr = np.asarray(data, dtype=float)
    if arr.ndim != 3:
        raise ValueError(f"Variable '{varname}' must be 3D; got shape {arr.shape}")

    shape = arr.shape
    lon_axes = [i for i, s in enumerate(shape) if s == nlon]
    lat_axes = [i for i, s in enumerate(shape) if s == nlat]
    time_axes = [i for i, s in enumerate(shape) if s == ntime]

    for iax in lon_axes:
        for jax in lat_axes:
            for kax in time_axes:
                if len({iax, jax, kax}) == 3:
                    return np.transpose(arr, (iax, jax, kax))

    raise ValueError(
        f"Unable to orient '{varname}' to (lon, lat, time). "
        f"shape={shape}, expected lengths lon={nlon}, lat={nlat}, time={ntime}"
    )


@dataclass(frozen=True)
class PercentileBinStats:
    bin_lower: np.ndarray
    bin_upper: np.ndarray
    bin_mass: np.ndarray
    bin_center: np.ndarray
    work_bin_avg: np.ndarray
    lift_bin_avg: np.ndarray
    ratio_bin: np.ndarray
    labels: List[str]


def build_percentile_bin_stats(
    percentile: np.ndarray,
    work_tail_avg: np.ndarray,
    lift_tail_avg: np.ndarray,
    work_total_avg: float,
    lift_total_avg: float,
) -> PercentileBinStats:
    p = np.asarray(percentile, dtype=float).reshape(-1)
    if p.size < 1:
        raise ValueError("No percentile thresholds found in the thresholded file")
    if np.any(np.diff(p) <= 0):
        raise ValueError("Percentile coordinate must be strictly increasing")
    if p[0] <= 0 or p[-1] >= 1:
        raise ValueError("Percentile coordinate values must be strictly inside (0, 1)")

    lower = np.concatenate(([0.0], p))
    upper = np.concatenate((p, [1.0]))
    mass = upper - lower
    nperc = p.size
    nbin = mass.size

    work_tail_avg = np.asarray(work_tail_avg, dtype=float).reshape(-1)
    lift_tail_avg = np.asarray(lift_tail_avg, dtype=float).reshape(-1)

    work_bin_avg = np.full(nbin, np.nan, dtype=float)
    lift_bin_avg = np.full(nbin, np.nan, dtype=float)

    tail_mass_first = 1.0 - p[0]
    work_bin_avg[0] = (work_total_avg - tail_mass_first * work_tail_avg[0]) / mass[0]
    lift_bin_avg[0] = (lift_total_avg - tail_mass_first * lift_tail_avg[0]) / mass[0]

    for ib in range(1, nperc):
        p_lo = p[ib - 1]
        p_hi = p[ib]
        tail_mass_lo = 1.0 - p_lo
        tail_mass_hi = 1.0 - p_hi
        work_bin_avg[ib] = (tail_mass_lo * work_tail_avg[ib - 1] - tail_mass_hi * work_tail_avg[ib]) / mass[ib]
        lift_bin_avg[ib] = (tail_mass_lo * lift_tail_avg[ib - 1] - tail_mass_hi * lift_tail_avg[ib]) / mass[ib]

    work_bin_avg[-1] = work_tail_avg[-1]
    lift_bin_avg[-1] = lift_tail_avg[-1]

    ratio_bin = lift_bin_avg / work_bin_avg
    ratio_bin[np.abs(work_bin_avg) < np.finfo(float).eps] = np.nan

    labels = [f"<P{_fmt(100.0 * p[0])}"]
    labels.extend(f"P{_fmt(100.0 * p[i - 1])}-P{_fmt(100.0 * p[i])}" for i in range(1, nperc))
    labels.append(f">P{_fmt(100.0 * p[-1])}")

    return PercentileBinStats(
        bin_lower=lower,
        bin_upper=upper,
        bin_mass=mass,
        bin_center=0.5 * (lower + upper),
        work_bin_avg=work_bin_avg,
        lift_bin_avg=lift_bin_avg,
        ratio_bin=ratio_bin,
        labels=labels,
    )


def _fmt(value: float) -> str:
    return (f"{value:.6f}").rstrip("0").rstrip(".")
