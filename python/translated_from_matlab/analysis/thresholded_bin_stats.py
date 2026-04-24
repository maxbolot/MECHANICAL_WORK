from dataclasses import dataclass

import numpy as np

from translated_from_matlab.lib.region_mask import build_region_mask
from translated_from_matlab.lib.thresholded_metrics import (
    PercentileBinStats,
    build_percentile_bin_stats,
    compute_thresholded_region_spatial_timeseries,
    compute_unthresholded_region_time_means,
)
from translated_from_matlab.lib.weights import weighted_nanmean


@dataclass(frozen=True)
class ThresholdedRunConfig:
    scenario: object
    region: object


@dataclass(frozen=True)
class ThresholdedBinStatsResult:
    stats: PercentileBinStats
    work_avg: np.ndarray
    lift_avg: np.ndarray


def compute_thresholded_bin_stats(cfg: ThresholdedRunConfig) -> ThresholdedBinStatsResult:
    ts = compute_thresholded_region_spatial_timeseries(cfg.scenario.thresholded_ncfile, build_region_mask_from_cfg(cfg))

    weights_result = cfg.scenario.time_weight_fn(ts.time, ts.time_units)
    time_weights_days = weights_result.weights_days

    if time_weights_days.size != ts.ntime:
        raise ValueError(
            f"Time-weight length ({time_weights_days.size}) does not match data time length ({ts.ntime})"
        )

    work_avg = np.full(ts.nperc, np.nan, dtype=float)
    lift_avg = np.full(ts.nperc, np.nan, dtype=float)
    for ip in range(ts.nperc):
        work_avg[ip] = weighted_nanmean(ts.work_spatial_avg[:, ip], time_weights_days)
        lift_avg[ip] = weighted_nanmean(ts.lift_spatial_avg[:, ip], time_weights_days)

    work_total_avg, lift_total_avg = compute_unthresholded_region_time_means(cfg, build_region_mask_from_cfg(cfg))

    stats = build_percentile_bin_stats(
        percentile=ts.percentile,
        work_tail_avg=work_avg,
        lift_tail_avg=lift_avg,
        work_total_avg=work_total_avg,
        lift_total_avg=lift_total_avg,
    )

    return ThresholdedBinStatsResult(stats=stats, work_avg=work_avg, lift_avg=lift_avg)


def build_region_mask_from_cfg(cfg):
    from translated_from_matlab.lib.io import read_thresholded_coordinates

    coords = read_thresholded_coordinates(cfg.scenario.thresholded_ncfile)
    return build_region_mask(coords.lat, coords.lon, cfg.region)
