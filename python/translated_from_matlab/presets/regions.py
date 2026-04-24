from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass(frozen=True)
class RegionConfig:
    name: str
    title_label: str
    result_label: str
    average_label: str
    lat_bounds: Optional[Tuple[float, float]]
    lon_bounds: Optional[Tuple[float, float]]


def region_tropics() -> RegionConfig:
    return RegionConfig(
        name="tropics",
        title_label="Tropical",
        result_label="Tropics (30°S to 30°N)",
        average_label="tropical band",
        lat_bounds=(-30.0, 30.0),
        lon_bounds=None,
    )
