# Project paths, shared by all R scripts.
# RAC_ROOT: directory holding data/ and results/ (default: the repository root).
# REPRO_GWAS_DIR: harmonised ANM/AMH/parental-lifespan GWAS and LD reference files.
# INSTRUMENT_DIR: metabolite GWAS summary statistics and extracted instruments (data/mqtl/).
RAC_ROOT <- Sys.getenv("RAC_ROOT", ".")
REPRO_GWAS_DIR <- Sys.getenv("REPRO_GWAS_DIR", file.path(RAC_ROOT, "external", "reproductive_gwas"))
INSTRUMENT_DIR <- Sys.getenv("INSTRUMENT_DIR", file.path(RAC_ROOT, "external", "metabolite_instruments"))
