from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import xarray as xr


@dataclass(frozen=True)
class ThresholdedCoordinates:
    lat: np.ndarray
    lon: np.ndarray
    time: np.ndarray
    percentile: np.ndarray
    time_units: Optional[str]


def ensure_input_exists(path: str) -> None:
    if not Path(path).is_file():
        raise FileNotFoundError(f"Input file not found: {path}")


def read_thresholded_coordinates(ncfile: str) -> ThresholdedCoordinates:
    ensure_input_exists(ncfile)
    with xr.open_dataset(ncfile, decode_times=False) as ds:
        lat = ds["lat"].values
        lon = ds["lon"].values
        time = ds["time"].values
        percentile = ds["percentile"].values
        time_units = ds["time"].attrs.get("units")
    return ThresholdedCoordinates(
        lat=np.asarray(lat, dtype=float),
        lon=np.asarray(lon, dtype=float),
        time=np.asarray(time, dtype=float),
        percentile=np.asarray(percentile, dtype=float),
        time_units=time_units,
    )


def read_standard_coordinates(ncfile: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Optional[str]]:
    ensure_input_exists(ncfile)
    with xr.open_dataset(ncfile, decode_times=False) as ds:
        lat = np.asarray(ds["lat"].values, dtype=float)
        lon = np.asarray(ds["lon"].values, dtype=float)
        time = np.asarray(ds["time"].values, dtype=float)
        time_units = ds["time"].attrs.get("units")
    return lat, lon, time, time_units


def read_thresholded_var(ncfile: str, varname: str) -> np.ndarray:
    ensure_input_exists(ncfile)
    with xr.open_dataset(ncfile, decode_times=False) as ds:
        if varname not in ds:
            raise KeyError(f"Variable '{varname}' not found in {ncfile}")
        arr = ds[varname].transpose("lon", "lat", "time", "percentile").values
    return np.asarray(arr, dtype=float)


def read_required_var(ncfile: str, varname: str) -> np.ndarray:
    ensure_input_exists(ncfile)
    with xr.open_dataset(ncfile, decode_times=False) as ds:
        if varname not in ds:
            raise KeyError(f"Variable '{varname}' not found in {ncfile}")
        da = ds[varname]
        data = _transpose_to_lon_lat_time(da).values
    return np.asarray(data, dtype=float)


def _transpose_to_lon_lat_time(da: xr.DataArray) -> xr.DataArray:
    """Return a view with canonical (lon, lat, time) order when possible."""
    dims = list(da.dims)

    lon_dim = _find_first_dim(dims, ("lon", "longitude", "grid_xt_coarse", "x"))
    lat_dim = _find_first_dim(dims, ("lat", "latitude", "grid_yt_coarse", "y"))
    time_dim = _find_first_dim(dims, ("time", "t"))

    if lon_dim is None or lat_dim is None or time_dim is None:
        raise ValueError(
            f"Variable '{da.name}' must include lon/lat/time dimensions; got {tuple(dims)}"
        )

    return da.transpose(lon_dim, lat_dim, time_dim)


def _find_first_dim(dims: List[str], candidates: Tuple[str, ...]) -> Optional[str]:
    for cand in candidates:
        if cand in dims:
            return cand
    return None
