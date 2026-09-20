"""Portable path configuration for PlantCellTwin workflows."""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("PCT_DATA_DIR", REPO_ROOT / "data")).expanduser().resolve()
RESULTS_DIR = Path(os.environ.get("PCT_RESULTS_DIR", REPO_ROOT / "results")).expanduser().resolve()

def require(path: Path, description: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {description}: {path}. See data/README.md or set PCT_DATA_DIR."
        )
    return path
