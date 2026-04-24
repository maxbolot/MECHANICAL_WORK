from dataclasses import dataclass
from typing import Callable
from typing import Optional

import numpy as np

from translated_from_matlab.lib.weights import (
    TimeWeightResult,
    compute_time_weights_control,
    compute_time_weights_plus4k,
)


@dataclass(frozen=True)
class ScenarioConfig:
    name: str
    label: str
    standard_ncfile: str
    standard_precip_ncfile: str
    thresholded_ncfile: str
    thresholded_precip_ncfile: str
    time_weight_fn: Callable[[np.ndarray, Optional[str]], TimeWeightResult]


def scenario_control() -> ScenarioConfig:
    return ScenarioConfig(
        name="control",
        label="Control",
        standard_ncfile="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180/work_2020010300_2022011200.nc",
        standard_precip_ncfile="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180/precip_2020010300_2022011200.nc",
        thresholded_ncfile="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_prate_thresholded/work_2020010300_2022011200.nc",
        thresholded_precip_ncfile="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_prate_thresholded/precip_2020010300_2022011200.nc",
        time_weight_fn=compute_time_weights_control,
    )


def scenario_plus4k() -> ScenarioConfig:
    return ScenarioConfig(
        name="plus4k",
        label="+4K",
        standard_ncfile="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv/work_2020010300_2022011800.nc",
        standard_precip_ncfile="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv/precip_2020010300_2022012000.nc",
        thresholded_ncfile="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded/work_2020010300_2022011800.nc",
        thresholded_precip_ncfile="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded/precip_2020010300_2022012000.nc",
        time_weight_fn=compute_time_weights_plus4k,
    )
