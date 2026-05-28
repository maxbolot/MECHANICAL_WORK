"""Path utilities for translated MATLAB workflows."""

from __future__ import annotations

import os
import platform
from pathlib import Path


def _looks_like_wsl() -> bool:
    if platform.system().lower() != "linux":
        return False

    release = platform.release().lower()
    if "microsoft" in release or "wsl" in release:
        return True

    # Fallback signal used in some WSL installs.
    return "WSL_DISTRO_NAME" in os.environ


def get_globalfv3_data_root() -> Path:
    """Resolve GLOBALFV3 root using env override, then platform defaults."""
    env_root = os.getenv("GLOBALFV3_DATA_ROOT")
    if env_root:
        return Path(env_root)

    if os.name == "nt":
        return Path(r"C:\climate_processed_data\GLOBALFV3")

    if _looks_like_wsl():
        return Path("/mnt/c/climate_processed_data/GLOBALFV3")

    return Path("/scratch/gpfs/mbolot/results/GLOBALFV3")
