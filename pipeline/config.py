"""Centralised configuration values.

Edit `PROJECT_DIR` to match your local checkout.
"""
from pathlib import Path

# Base project directory
PROJECT_DIR: Path = Path.home() / "git" / "valleyfevermutation"

# I/O defaults
DATA_DIR: Path = PROJECT_DIR / "dataset" / "proteins"
MD_SIM_DIR: Path = PROJECT_DIR / "ai"

# Model / graph hyper‑parameters
ESM_MODEL_NAME: str = "esm2_t33_650M_UR50D"
ESM_LAYER: int = 33
KNN_K: int = 16
DISTANCE_CUTOFF: tuple[float, float] = (4.0, 8.0)  # Å
