from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from translated_from_matlab.lib.io import ensure_input_exists, read_required_var, read_standard_coordinates
from translated_from_matlab.lib.weights import compute_time_weights_control, compute_time_weights_plus4k, weighted_nanmean


@dataclass(frozen=True)
class WorkLiftRegion:
    name: str
    latmin: float
    latmax: float
    lonmin: Optional[float]
    lonmax: Optional[float]


@dataclass(frozen=True)
class WorkLiftAllRegionsConfig:
    control_file: str
    warming_file: str
    regions: Sequence[WorkLiftRegion]
    make_plots: bool
    print_summary: bool


@dataclass(frozen=True)
class RegionSummary:
    work: np.ndarray
    lift: np.ndarray
    ke: np.ndarray
    ratio: np.ndarray
    missing_steps: np.ndarray


@dataclass(frozen=True)
class WorkLiftAllRegionsResult:
    control: RegionSummary
    warming: RegionSummary
    newY: np.ndarray
    regions: Sequence[WorkLiftRegion]
    figure: Optional[plt.Figure]


def default_regions() -> List[WorkLiftRegion]:
    return [
        WorkLiftRegion(name="global", latmin=-90.0, latmax=90.0, lonmin=None, lonmax=None),
        WorkLiftRegion(name="tropics", latmin=-30.0, latmax=30.0, lonmin=None, lonmax=None),
        WorkLiftRegion(name="maritime continent", latmin=-15.0, latmax=15.0, lonmin=90.0, lonmax=150.0),
        WorkLiftRegion(name="north. midlatitudes", latmin=30.0, latmax=60.0, lonmin=None, lonmax=None),
    ]


def apply_defaults(cfg: Optional[dict] = None) -> WorkLiftAllRegionsConfig:
    if cfg is None:
        cfg = {}

    control_file = cfg.get(
        "control_file",
        "/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180/work_2020010300_2022011200.nc",
    )
    warming_file = cfg.get(
        "warming_file",
        "/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv/work_2020010300_2022011800.nc",
    )
    regions = cfg.get("regions", default_regions())
    make_plots = bool(cfg.get("make_plots", True))
    print_summary = bool(cfg.get("print_summary", True))

    return WorkLiftAllRegionsConfig(
        control_file=control_file,
        warming_file=warming_file,
        regions=regions,
        make_plots=make_plots,
        print_summary=print_summary,
    )


def run_work_lift_all_regions_analysis(cfg: Optional[dict] = None) -> WorkLiftAllRegionsResult:
    cfg_obj = apply_defaults(cfg)

    ensure_input_exists(cfg_obj.control_file)
    ensure_input_exists(cfg_obj.warming_file)

    print("Reading control simulation:", cfg_obj.control_file)
    control = process_work_lift_regions_file(cfg_obj.control_file, cfg_obj.regions, "control")

    print("Reading warming simulation:", cfg_obj.warming_file)
    warming = process_work_lift_regions_file(cfg_obj.warming_file, cfg_obj.regions, "warming")

    y = np.zeros((8, 2), dtype=float)
    for i in range(4):
        y[2 * i, :] = [control.lift[i], control.ke[i]]
        y[2 * i + 1, :] = [warming.lift[i], warming.ke[i]]

    newY = np.zeros((12, 2), dtype=float)
    newY[[0, 1, 3, 4, 6, 7, 9, 10], :] = y

    fig = None
    if cfg_obj.make_plots:
        fig = plot_legacy_bars(newY)

    if cfg_obj.print_summary:
        print("\n=== Spatial + Time Averaged Results ===")
        print(f"{'Region':<22} {'Simulation':<12} {'work':>12} {'lift':>12} {'ke':>12} {'lift/work':>12} {'missing':>10}")
        for i, region in enumerate(cfg_obj.regions):
            print(
                f"{region.name:<22} {'control':<12} {control.work[i]:12.6f} {control.lift[i]:12.6f} {control.ke[i]:12.6f} {control.ratio[i]:12.6f} {int(control.missing_steps[i]):10d}"
            )
            print(
                f"{region.name:<22} {'warming':<12} {warming.work[i]:12.6f} {warming.lift[i]:12.6f} {warming.ke[i]:12.6f} {warming.ratio[i]:12.6f} {int(warming.missing_steps[i]):10d}"
            )

    return WorkLiftAllRegionsResult(
        control=control,
        warming=warming,
        newY=newY,
        regions=cfg_obj.regions,
        figure=fig,
    )


def process_work_lift_regions_file(ncfile: str, regions: Sequence[WorkLiftRegion], simulation_name: str) -> RegionSummary:
    lat, lon, time, time_units = read_standard_coordinates(ncfile)
    work3d = read_required_var(ncfile, "work")
    lift3d = read_required_var(ncfile, "lift")

    if simulation_name.lower() in ("warming", "plus4k"):
        tw = compute_time_weights_plus4k(time, time_units)
    elif simulation_name.lower() == "control":
        tw = compute_time_weights_control(time, time_units)
    else:
        raise ValueError("Unknown simulation name: %s" % simulation_name)

    time_weights_days = tw.weights_days
    missing_steps = tw.missing_steps

    nreg = len(regions)
    out_work = np.zeros(nreg, dtype=float)
    out_lift = np.zeros(nreg, dtype=float)
    out_ke = np.zeros(nreg, dtype=float)
    out_ratio = np.zeros(nreg, dtype=float)
    out_missing = np.zeros(nreg, dtype=float)

    lon_wrapped = np.mod(lon.reshape(-1), 360.0)

    for ir, r in enumerate(regions):
        lat_idx = np.where((lat >= r.latmin) & (lat <= r.latmax))[0]
        if r.lonmin is None or r.lonmax is None:
            lon_idx = np.arange(lon.size)
        else:
            lon_idx = np.where((lon_wrapped >= r.lonmin) & (lon_wrapped <= r.lonmax))[0]

        if lat_idx.size == 0 or lon_idx.size == 0:
            raise ValueError("No grid cells found for region: %s" % r.name)

        lat_weights = np.cos(np.deg2rad(lat[lat_idx]))
        w3 = lat_weights.reshape(1, -1, 1)

        wr = work3d[np.ix_(lon_idx, lat_idx, np.arange(work3d.shape[2]))]
        lr = lift3d[np.ix_(lon_idx, lat_idx, np.arange(lift3d.shape[2]))]

        w_num = np.nansum(wr * w3, axis=(0, 1))
        w_den = np.nansum((~np.isnan(wr)).astype(float) * w3, axis=(0, 1))
        work_ts = w_num / w_den

        l_num = np.nansum(lr * w3, axis=(0, 1))
        l_den = np.nansum((~np.isnan(lr)).astype(float) * w3, axis=(0, 1))
        lift_ts = l_num / l_den

        work_avg = weighted_nanmean(work_ts, time_weights_days)
        lift_avg = weighted_nanmean(lift_ts, time_weights_days)

        out_work[ir] = work_avg
        out_lift[ir] = lift_avg
        out_ke[ir] = work_avg - lift_avg
        out_ratio[ir] = np.nan if abs(work_avg) < np.finfo(float).eps else lift_avg / work_avg
        out_missing[ir] = missing_steps

    return RegionSummary(work=out_work, lift=out_lift, ke=out_ke, ratio=out_ratio, missing_steps=out_missing)


def plot_legacy_bars(newY: np.ndarray) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8, 3.8), dpi=110)

    x = np.arange(newY.shape[0])
    lift = newY[:, 0]
    ke = newY[:, 1]

    control_dark = (0.10, 0.35, 0.80)
    control_light = (0.65, 0.80, 0.98)
    warming_dark = (0.85, 0.20, 0.20)
    warming_light = (0.98, 0.72, 0.72)
    white_bar = (1.0, 1.0, 1.0)

    c1 = np.array([
        control_dark, warming_dark, white_bar, control_dark, warming_dark, white_bar,
        control_dark, warming_dark, white_bar, control_dark, warming_dark, white_bar,
    ])
    c2 = np.array([
        control_light, warming_light, white_bar, control_light, warming_light, white_bar,
        control_light, warming_light, white_bar, control_light, warming_light, white_bar,
    ])

    b1 = ax.bar(x, lift, color=c1, width=0.8, label="lift")
    b2 = ax.bar(x, ke, bottom=lift, color=c2, width=0.8, label="ke = work - lift")

    ax.set_xlim(-0.5, 11.5)
    ax.set_xticks([1.5, 4.5, 7.5, 10.5])
    ax.set_xticklabels(["global", "tropics", "maritime continent", "north. midlatitudes"])
    ax.set_ylabel("W m$^{-2}$")
    ax.set_title("Legacy Axis 1: work = lift + ke")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper left", frameon=False)

    fig.tight_layout()
    return fig


def build_summary_table(result: WorkLiftAllRegionsResult) -> pd.DataFrame:
    rows = []
    for i, region in enumerate(result.regions):
        rows.append(
            {
                "region": region.name,
                "simulation": "control",
                "work": result.control.work[i],
                "lift": result.control.lift[i],
                "ke": result.control.ke[i],
                "lift_over_work": result.control.ratio[i],
                "missing_steps": int(result.control.missing_steps[i]),
            }
        )
        rows.append(
            {
                "region": region.name,
                "simulation": "warming",
                "work": result.warming.work[i],
                "lift": result.warming.lift[i],
                "ke": result.warming.ke[i],
                "lift_over_work": result.warming.ratio[i],
                "missing_steps": int(result.warming.missing_steps[i]),
            }
        )
    return pd.DataFrame(rows)
