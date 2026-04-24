from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import xarray as xr

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
    HAS_CARTOPY = True
except Exception:
    HAS_CARTOPY = False

try:
    from mpl_toolkits.basemap import Basemap
    HAS_BASEMAP = True
except Exception:
    HAS_BASEMAP = False

from translated_from_matlab.lib.io import ensure_input_exists
from translated_from_matlab.lib.weights import normalize_time_axis_to_days
from translated_from_matlab.presets.scenarios import scenario_control, scenario_plus4k


@dataclass(frozen=True)
class WorkRatioLonLatConfig:
    simulation: str
    temporal_mode: str
    instantaneous_date: datetime
    range_start: datetime
    range_end: datetime
    threshold_percentile: float
    make_plots: bool


@dataclass(frozen=True)
class WorkRatioLonLatResult:
    ncfile: str
    simulation: str
    simulation_label: str
    is_thresholded: bool
    temporal_mode: str
    time_indices: np.ndarray
    time_label: str
    lon: np.ndarray
    lat: np.ndarray
    work_map: np.ndarray
    lift_map: np.ndarray
    ratio_map: np.ndarray
    selected_percentile: Optional[float]
    figure: Optional[plt.Figure]


def apply_defaults(cfg: Optional[dict] = None) -> WorkRatioLonLatConfig:
    if cfg is None:
        cfg = {}

    return WorkRatioLonLatConfig(
        simulation=cfg.get("simulation", "control"),
        temporal_mode=cfg.get("temporal_mode", "instantaneous"),
        instantaneous_date=cfg.get("instantaneous_date", datetime(2021, 1, 1, 0, 0, 0, tzinfo=timezone.utc)),
        range_start=cfg.get("range_start", datetime(2020, 8, 1, 0, 0, 0, tzinfo=timezone.utc)),
        range_end=cfg.get("range_end", datetime(2020, 8, 30, 23, 59, 59, tzinfo=timezone.utc)),
        threshold_percentile=float(cfg.get("threshold_percentile", 0.5)),
        make_plots=bool(cfg.get("make_plots", True)),
    )


def run_work_ratio_lonlat_map_analysis(cfg: Optional[dict] = None) -> WorkRatioLonLatResult:
    cfg_obj = apply_defaults(cfg)

    ncfile, simulation_label, is_thresholded = resolve_work_ratio_input_file(cfg_obj.simulation)
    lon, lat, time_axis, time_units = read_coordinates(ncfile)
    time_dt = netcdf_time_to_datetime(time_axis, time_units)

    time_indices, time_label = select_time_indices(
        time_dt,
        cfg_obj.temporal_mode,
        cfg_obj.instantaneous_date,
        cfg_obj.range_start,
        cfg_obj.range_end,
    )

    if is_thresholded:
        work_data, selected_percentile = read_lon_lat_time_from_var(ncfile, "work", cfg_obj.threshold_percentile)
        lift_data, _ = read_lon_lat_time_from_var(ncfile, "lift", cfg_obj.threshold_percentile)
        simulation_label = "%s (thresholded p%g)" % (simulation_label, 100.0 * selected_percentile)
    else:
        work_data, _ = read_lon_lat_time_from_var(ncfile, "work", np.nan)
        lift_data, _ = read_lon_lat_time_from_var(ncfile, "lift", np.nan)
        selected_percentile = None

    work_map = mean_over_time(work_data, time_indices)
    lift_map = mean_over_time(lift_data, time_indices)

    ratio_map = lift_map / work_map
    ratio_map[np.abs(work_map) < np.finfo(float).eps] = np.nan

    fig = None
    if cfg_obj.make_plots:
        fig = plot_maps(lon, lat, work_map, ratio_map, simulation_label, time_label)

    return WorkRatioLonLatResult(
        ncfile=ncfile,
        simulation=cfg_obj.simulation,
        simulation_label=simulation_label,
        is_thresholded=is_thresholded,
        temporal_mode=cfg_obj.temporal_mode,
        time_indices=time_indices,
        time_label=time_label,
        lon=lon,
        lat=lat,
        work_map=work_map,
        lift_map=lift_map,
        ratio_map=ratio_map,
        selected_percentile=selected_percentile,
        figure=fig,
    )


def resolve_work_ratio_input_file(simulation: str) -> Tuple[str, str, bool]:
    sim = simulation.lower()
    if sim == "control":
        sc = scenario_control()
        ncfile = sc.standard_ncfile
        label = "Control"
        is_thresholded = False
    elif sim == "warming":
        sc = scenario_plus4k()
        ncfile = sc.standard_ncfile
        label = "Warming (+4K)"
        is_thresholded = False
    elif sim == "control_prate_thresholded":
        sc = scenario_control()
        ncfile = sc.thresholded_ncfile
        label = "Control"
        is_thresholded = True
    elif sim == "warming_prate_thresholded":
        sc = scenario_plus4k()
        ncfile = sc.thresholded_ncfile
        label = "Warming (+4K)"
        is_thresholded = True
    else:
        raise ValueError("Unsupported simulation: %s" % simulation)

    ensure_input_exists(ncfile)
    return ncfile, label, is_thresholded


def read_coordinates(ncfile: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Optional[str]]:
    with xr.open_dataset(ncfile, decode_times=False) as ds:
        lon = np.asarray(ds["lon"].values, dtype=float)
        lat = np.asarray(ds["lat"].values, dtype=float)
        time_axis = np.asarray(ds["time"].values, dtype=float)
        time_units = ds["time"].attrs.get("units")
    return lon, lat, time_axis, time_units


def read_lon_lat_time_from_var(ncfile: str, varname: str, target_percentile: float) -> Tuple[np.ndarray, Optional[float]]:
    with xr.open_dataset(ncfile, decode_times=False) as ds:
        if varname not in ds:
            raise KeyError("Variable '%s' not found in %s" % (varname, ncfile))

        da = ds[varname]
        if "percentile" in da.dims:
            pvals = np.asarray(ds["percentile"].values, dtype=float)
            if np.isfinite(target_percentile):
                pidx = int(np.argmin(np.abs(pvals - float(target_percentile))))
            else:
                pidx = 0
            selected_percentile = float(pvals[pidx])
            da = da.transpose("lon", "lat", "time", "percentile").isel(percentile=pidx)
            return np.asarray(da.values, dtype=float), selected_percentile

        da = da.transpose("lon", "lat", "time")
        return np.asarray(da.values, dtype=float), None


def netcdf_time_to_datetime(time_axis: np.ndarray, time_units: Optional[str]) -> pd.DatetimeIndex:
    time_days = normalize_time_axis_to_days(time_axis, time_units)
    seconds = np.asarray(time_days, dtype=float) * 86400.0
    return pd.to_datetime(seconds, unit="s", utc=True)


def select_time_indices(
    time_dt: pd.DatetimeIndex,
    temporal_mode: str,
    inst_date: datetime,
    range_start: datetime,
    range_end: datetime,
) -> Tuple[np.ndarray, str]:
    mode = temporal_mode.lower()

    if inst_date.tzinfo is None:
        inst_date = inst_date.replace(tzinfo=timezone.utc)
    if range_start.tzinfo is None:
        range_start = range_start.replace(tzinfo=timezone.utc)
    if range_end.tzinfo is None:
        range_end = range_end.replace(tzinfo=timezone.utc)

    if mode == "instantaneous":
        target = pd.Timestamp(inst_date)
        idx = int(np.argmin(np.abs((time_dt - target).total_seconds())))
        label = "instantaneous %s" % time_dt[idx].strftime("%Y-%m-%d %H:%M:%S UTC")
        return np.asarray([idx], dtype=int), label

    if mode == "timerange":
        mask = (time_dt >= pd.Timestamp(range_start)) & (time_dt <= pd.Timestamp(range_end))
        idx = np.where(mask)[0]
        if idx.size == 0:
            raise ValueError("No timestamps found inside selected range")
        label = "time-mean %s to %s" % (time_dt[idx[0]].strftime("%Y-%m-%d"), time_dt[idx[-1]].strftime("%Y-%m-%d"))
        return idx, label

    raise ValueError("Unsupported temporal_mode: %s" % temporal_mode)


def mean_over_time(field_data: np.ndarray, time_indices: Sequence[int]) -> np.ndarray:
    arr = np.asarray(field_data, dtype=float)
    if arr.ndim != 3:
        raise ValueError("Expected field_data dimensions [lon, lat, time], got %s" % (arr.shape,))
    return np.nanmean(arr[:, :, np.asarray(time_indices, dtype=int)], axis=2)


def plot_maps(
    lon: np.ndarray,
    lat: np.ndarray,
    work_map: np.ndarray,
    ratio_map: np.ndarray,
    simulation_label: str,
    time_label: str,
) -> plt.Figure:
    lon_plot, work_plot, ratio_plot = wrap_longitude_for_plot(lon, work_map, ratio_map)
    x_grid, y_grid = np.meshgrid(lon_plot.astype(float), lat.astype(float), indexing="ij")

    if HAS_CARTOPY:
        fig = plt.figure(figsize=(13, 5), dpi=110, constrained_layout=True)
        ax1 = fig.add_subplot(2, 1, 1, projection=ccrs.PlateCarree())
        ax2 = fig.add_subplot(2, 1, 2, projection=ccrs.PlateCarree())
        axes = [ax1, ax2]
        transform = ccrs.PlateCarree()
    else:
        fig, axes = plt.subplots(2, 1, figsize=(13, 5), dpi=110, constrained_layout=True)
        transform = None

    if transform is None:
        im1 = axes[0].pcolormesh(x_grid, y_grid, work_plot, shading="auto", cmap="viridis", vmin=0, vmax=50)
    else:
        im1 = axes[0].pcolormesh(
            x_grid,
            y_grid,
            work_plot,
            shading="auto",
            cmap="viridis",
            vmin=0,
            vmax=50,
            transform=transform,
        )

    _draw_coastlines(axes[0])
    _draw_map_coordinates(axes[0])
    axes[0].set_title("a. Work (%s)" % time_label)
    axes[0].set_xlim(-180, 180)
    axes[0].set_ylim(-60, 60)
    axes[0].set_xlabel("Longitude")
    axes[0].set_ylabel("Latitude")
    fig.colorbar(im1, ax=axes[0], location="right", label="W m$^{-2}$")

    if transform is None:
        im2 = axes[1].pcolormesh(x_grid, y_grid, ratio_plot, shading="auto", cmap="viridis", vmin=0, vmax=1)
    else:
        im2 = axes[1].pcolormesh(
            x_grid,
            y_grid,
            ratio_plot,
            shading="auto",
            cmap="viridis",
            vmin=0,
            vmax=1,
            transform=transform,
        )

    _draw_coastlines(axes[1])
    _draw_map_coordinates(axes[1])
    axes[1].set_title("b. Lift/Work (%s)" % time_label)
    axes[1].set_xlim(-180, 180)
    axes[1].set_ylim(-60, 60)
    axes[1].set_xlabel("Longitude")
    axes[1].set_ylabel("Latitude")
    fig.colorbar(im2, ax=axes[1], location="right", label="unitless")

    fig.suptitle("%s: Work and Lift/Work in Lon-Lat Coordinates" % simulation_label)
    return fig


def wrap_longitude_for_plot(lon: np.ndarray, work_map: np.ndarray, ratio_map: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Wrap longitudes from [0, 360) style grids into [-180, 180] for plotting."""
    lon_arr = np.asarray(lon, dtype=float).reshape(-1)
    work_arr = np.asarray(work_map, dtype=float)
    ratio_arr = np.asarray(ratio_map, dtype=float)

    if work_arr.shape[0] != lon_arr.size or ratio_arr.shape[0] != lon_arr.size:
        raise ValueError(
            "Longitude length does not match map first dimension: "
            "len(lon)=%d, work_shape=%s, ratio_shape=%s"
            % (lon_arr.size, work_arr.shape, ratio_arr.shape)
        )

    lon_wrapped = ((lon_arr + 180.0) % 360.0) - 180.0
    sort_idx = np.argsort(lon_wrapped)

    return lon_wrapped[sort_idx], work_arr[sort_idx, :], ratio_arr[sort_idx, :]


def _draw_coastlines(ax) -> None:
    """Draw coastlines if cartopy or basemap is available."""
    if HAS_CARTOPY:
        ax.coastlines(resolution="110m", linewidth=0.5, color="black")
        try:
            ax.add_feature(cfeature.BORDERS, linewidth=0.25, edgecolor="black")
        except Exception:
            pass
        return

    if HAS_BASEMAP:
        m = Basemap(
            projection="cyl",
            llcrnrlon=-180,
            urcrnrlon=180,
            llcrnrlat=-60,
            urcrnrlat=60,
            ax=ax,
            resolution="c",
        )
        m.drawcoastlines(linewidth=0.5, color="black")


def _draw_map_coordinates(ax) -> None:
    """Draw labeled lon/lat coordinate lines for map readability."""
    xticks = np.arange(-180, 181, 60)
    yticks = np.arange(-60, 61, 20)

    if HAS_CARTOPY:
        gridliner = ax.gridlines(
            crs=ccrs.PlateCarree(),
            draw_labels=True,
            linewidth=0.4,
            color="gray",
            alpha=0.5,
            linestyle="--",
        )
        gridliner.top_labels = False
        gridliner.right_labels = False
        gridliner.xlocator = mticker.FixedLocator(xticks)
        gridliner.ylocator = mticker.FixedLocator(yticks)
        gridliner.xformatter = LONGITUDE_FORMATTER
        gridliner.yformatter = LATITUDE_FORMATTER
        gridliner.xlabel_style = {"size": 9}
        gridliner.ylabel_style = {"size": 9}
        return

    if HAS_BASEMAP:
        m = Basemap(
            projection="cyl",
            llcrnrlon=-180,
            urcrnrlon=180,
            llcrnrlat=-60,
            urcrnrlat=60,
            ax=ax,
            resolution="c",
        )
        m.drawparallels(yticks, labels=[1, 0, 0, 0], linewidth=0.4, dashes=[2, 2], color="gray", fontsize=9)
        m.drawmeridians(xticks, labels=[0, 0, 0, 1], linewidth=0.4, dashes=[2, 2], color="gray", fontsize=9)
        return

    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.grid(True, which="major", linestyle="--", linewidth=0.4, alpha=0.5)
