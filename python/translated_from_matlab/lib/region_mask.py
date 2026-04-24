from dataclasses import dataclass

import numpy as np

from translated_from_matlab.presets.regions import RegionConfig


@dataclass(frozen=True)
class RegionMask:
    lat_idx: np.ndarray
    lon_idx: np.ndarray
    lat_weights: np.ndarray
    selection_message: str


def build_region_mask(lat: np.ndarray, lon: np.ndarray, region: RegionConfig) -> RegionMask:
    lat_values = np.asarray(lat, dtype=float).reshape(-1)
    lon_values = np.asarray(lon, dtype=float).reshape(-1)
    lon_wrapped = np.mod(lon_values, 360.0)

    if region.lat_bounds is None:
        lat_idx = np.arange(lat_values.size)
    else:
        lat_idx = np.where((lat_values >= region.lat_bounds[0]) & (lat_values <= region.lat_bounds[1]))[0]

    if region.lon_bounds is None:
        lon_idx = np.arange(lon_values.size)
    else:
        lon_idx = np.where((lon_wrapped >= region.lon_bounds[0]) & (lon_wrapped <= region.lon_bounds[1]))[0]

    if lat_idx.size == 0 or lon_idx.size == 0:
        raise ValueError(f"No grid points found for region {region.name}")

    lat_selected = lat_values[lat_idx]
    lat_weights = np.cos(np.deg2rad(lat_selected))
    if np.all(lat_weights == 0):
        raise ValueError(f"Latitude weights are all zero in the selected {region.name} region")

    if region.name == "tropics":
        msg = f"Tropical latitude band selected: {lat_selected.min()} to {lat_selected.max()} degrees"
    else:
        msg = f"Region selected: {region.name}"

    return RegionMask(
        lat_idx=lat_idx,
        lon_idx=lon_idx,
        lat_weights=lat_weights.astype(float),
        selection_message=msg,
    )
