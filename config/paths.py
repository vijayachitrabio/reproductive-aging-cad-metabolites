"""Project paths, shared by all Python scripts (see config/paths.R for the meaning of each)."""
import os
from pathlib import Path

RAC_ROOT = Path(os.environ.get("RAC_ROOT", "."))
REPRO_GWAS_DIR = Path(os.environ.get("REPRO_GWAS_DIR", RAC_ROOT / "external" / "reproductive_gwas"))
INSTRUMENT_DIR = Path(os.environ.get("INSTRUMENT_DIR", RAC_ROOT / "external" / "metabolite_instruments"))
