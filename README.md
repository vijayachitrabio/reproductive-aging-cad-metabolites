# Shared genetic and metabolic architecture of reproductive aging, coronary artery disease, and parental lifespan

Analysis code for the article:

> Modhukur V, Khatun M, Patel N, Destouni A, Lingasamy P, Gromiha MM, Salumets A. *Shared genetic and metabolic architecture of reproductive aging, coronary artery disease, and parental lifespan: a Mendelian randomization study.* Manuscript, 2026.

## Overview

The study asks whether circulating metabolic traits link age at natural menopause (ANM) with coronary artery disease (CAD) and parental lifespan, and whether these relationships reflect causal mediation or shared genetic architecture. The analyses are:

1. Two-sample Mendelian randomization (MR) screens of 56 NMR metabolites against ANM and CAD, with FDR correction, and follow-up of prioritized metabolites against anti-Müllerian hormone (AMH) and parental lifespan
2. Reverse MR, pleiotropy-robust estimators (MR-PRESSO, MR-RAPS, contamination mixture, leave-one-out), Steiger directionality and LD-clumping sensitivity analyses
3. Multivariable MR of a four-trait metabolic axis
4. Total-effect MR of ANM on CAD and on parental lifespan
5. Genome-wide genetic correlation (LDSC)
6. Colocalization (coloc.abf) and multi-signal fine-mapping (SuSiE)
7. Expression and regulatory analyses at the FADS1–FADS2 locus (S-PrediXcan TWAS, GTEx v8 cis-eQTL lookups, AlphaGenome predictions)
8. STRING protein-interaction network and enrichment analyses
9. Exploratory FADS1 docking (AutoDock Vina) and 100-ns molecular dynamics with MM-GBSA (GROMACS, gmx_MMPBSA)

The study uses only publicly available GWAS summary statistics and public reference resources. It does not use individual-level participant data. No data files are redistributed in this repository; `data/`, `results/` and `external/` are excluded from version control.

## Data sources

| Data | Source |
|---|---|
| 56 NMR metabolites (UK Biobank, 246,683 participants) | Alasoo K, Rahu I. Zenodo, 2024. <https://doi.org/10.5281/zenodo.13821209> |
| Age at natural menopause | Ruth KS et al. Nature 2021;596:393–397. <https://doi.org/10.1038/s41586-021-03779-7> |
| Coronary artery disease | Nikpay M et al. Nat Genet 2015;47:1121–1130 (OpenGWAS `ieu-a-7`). <https://doi.org/10.1038/ng.3396> |
| Parental lifespan | Timmers PRHJ et al. eLife 2019;8:e39856. <https://doi.org/10.7554/eLife.39856> |
| Anti-Müllerian hormone | Pujol-Gualdo N et al. Hum Reprod 2024;39:1564–1572 (GWAS Catalog `GCST90428625`). <https://doi.org/10.1093/humrep/deae117> |
| Gene expression | GTEx v8 (<https://gtexportal.org>); MASHR S-PrediXcan models |
| FADS1 structure | AlphaFold Protein Structure Database model for UniProt O60427 (<https://www.uniprot.org/uniprotkb/O60427/entry>) |
| Ligands | PubChem CIDs 135539077 (Luminespib), 11327430 (NU7441), 59416313 (D5D-IN-326) |

The LDSC, TWAS and SuSiE analyses read dense summary statistics for CAD (`ieu-a-7`) and the prioritized metabolites (OpenGWAS `met-d-*` datasets) from `data/dense_sumstats/`.

## Repository layout

```
config/            paths.R, paths.py: project paths read from environment variables
environment/       software and package versions
scripts/
  01_instruments/                 metabolite instrument extraction
  02_forward_mr/                  metabolite -> ANM, CAD, AMH, parental lifespan; modules and enrichment
  03_reverse_mr/                  ANM, CAD, parental lifespan -> prioritized metabolites
  04_sensitivity/                 MR-PRESSO, MR-RAPS/ConMix/leave-one-out, Steiger, LD clumping
  05_mvmr/                        multivariable MR of the metabolic axis
  06_total_effect_mr/             ANM -> CAD (primary and expanded instruments), ANM -> parental lifespan
  07_ldsc/                        formatting, munging and LDSC genetic correlations
  08_colocalization_finemapping/  coloc.abf and SuSiE colocalization
  09_expression_regulatory/       S-PrediXcan TWAS, GTEx eQTL lookups, AlphaGenome
  10_network_enrichment/          STRING network, robustness, restricted-background enrichment
  11_evidence_integration/        evidence integration and tiers
  12_docking/                     AutoSite pocket, Vina configurations, receptor and ligand inputs
  13_molecular_dynamics/          GROMACS parameter files, system setup, post-MD analysis, MM-GBSA
  14_figures/                     Figures 2–5 and 7, Supplementary Figures S1–S7
```

## Configuration

Run all scripts from the repository root. Paths are set with environment variables (defaults in brackets):

| Variable | Meaning |
|---|---|
| `RAC_ROOT` [`.`] | directory holding `data/` and `results/` |
| `REPRO_GWAS_DIR` [`external/reproductive_gwas`] | harmonised ANM, AMH and parental-lifespan summary statistics, ANM lead loci, LDSC reference files (`data/ldsc_ref/`) and LDSC munging tool (`tools/ldsc/`) |
| `INSTRUMENT_DIR` [`external/metabolite_instruments`] | metabolite GWAS (`data/mqtl/mqtl_full_gwas/`) and extracted instruments (`data/mqtl/mqtl_gwas_instruments.csv`) |
| `TWAS_TOOL_DIR` [`external/metaxcan`] | MetaXcan, MASHR models and SNP covariance file |
| `LDSC_DIR` [`external/ldsc`] | LDSC software |
| `OPENGWAS_JWT` | OpenGWAS token, needed by scripts that query OpenGWAS |
| `ALPHA_GENOME_API_KEY` | AlphaGenome API key |
| `MD_ROOT` [`$HOME`] | parent directory of the per-system MD run folders |

## Pipeline

| Step | Script | Manuscript output |
|---|---|---|
| Instruments | `01_instruments/extract_metabolite_instruments.R` | instruments for all MR analyses |
| Metabolite → ANM | `02_forward_mr/metabolites_to_anm_mr.R` | Figure 2A, Table S19 |
| Metabolite → CAD | `02_forward_mr/metabolites_to_cad_mr_all.R`, `metabolites_to_cad_mr_targeted.R`, `summarize_shared_metabolites.R` | Figure 2B, Table S19 |
| Metabolite → AMH | `02_forward_mr/metabolites_to_amh_mr.R` | Table S19, Figure S5 |
| Metabolite → parental lifespan | `02_forward_mr/metabolites_to_lifespan_mr.R` | Table S9 |
| Modules and enrichment | `02_forward_mr/assign_metabolite_modules.R`, `module_enrichment_test.R` | Table 1, Table S2 |
| Reverse MR | `03_reverse_mr/*.R` | Table S20, Table S9, Figure S1 |
| Sensitivity | `04_sensitivity/ld_clumping_sensitivity.R`, `mr_presso.R`, `mr_raps_conmix_leave_one_out.R`, `steiger_directionality.R` | Tables S1, S3a, S3b |
| MVMR | `05_mvmr/metabolic_axis_mvmr.R` | Figure 3, Table S21 |
| Total-effect MR | `06_total_effect_mr/*.R` | Tables S22a–c, Table S9 |
| LDSC | `07_ldsc/format_and_munge_sumstats.py`, `run_ldsc_rg.sh`, `summarize_ldsc_rg.R` | Figure 4, Table 3, Table S5 |
| Colocalization | `08_colocalization_finemapping/metabolite_ovarian_coloc.R`, `metabolite_cad_coloc.R` | Table S23, Figures S3–S4 |
| Fine-mapping | `08_colocalization_finemapping/susie_coloc_cad_loci.R`, `susie_coloc_fads_apoe.R` | Table 2, Tables S4a–c |
| TWAS | `09_expression_regulatory/prepare_twas_sumstats.py`, `run_twas.sh`, `compile_twas_results.py` | Table S8 |
| GTEx eQTL | `09_expression_regulatory/query_gtex_eqtl.py` | Table S6, Table 6 |
| AlphaGenome | `09_expression_regulatory/run_alphagenome.sh`, `alphagenome_summary_table.py` | Tables S7a–c, Table 6 |
| Network | `10_network_enrichment/string_ppi_network.R`, `ppi_robustness_overlay.R`, `restricted_background_enrichment.py` | Figure 5, Table 4, Tables S10–S15 |
| Evidence integration | `11_evidence_integration/evidence_tiers.R`, `evidence_integration.R` | Figure S2 |
| Docking | `12_docking/` | Figure 6, Table 5, Table S16B |
| Molecular dynamics | `13_molecular_dynamics/` | Figure 7, Table 5, Tables S16A, S17, S18, Figures S6–S7 |
| Figures | `14_figures/main_figures_2_and_5.R`, `figure3_mvmr.R`, `figure4_ldsc.R`, `supplementary_figures_S1_S5.R` | Figures 2–5, S1–S5 |
| MD figures | `14_figures/md_figures_7_S6_S7.py` (script by Naisarg Patel) | Figure 7, Figures S6–S7 |

Figure 1 (study design) was drawn in BioRender. Figure 6 was drawn from the docking poses. Figures 7, S6 and S7 are plotted from the GROMACS outputs by `14_figures/md_figures_7_S6_S7.py`, written by Naisarg Patel.

## Docking

`12_docking/` contains the AutoSite pocket summary, the AutoDock Vina configuration for each ligand (search box centre (−1.228518, 8.457491, −13.414805) Å, 21 × 21 × 21 Å, exhaustiveness 32), the AlphaFold FADS1 model (`inputs/FADS1.pdb`), the receptor PDBQT prepared with ADFR suite `prepare_receptor4`, and the ligand PDBQT files prepared with Meeko 0.7.1 from the PubChem structures (C1 = Luminespib, C2 = NU7441, T = D5D-IN-326). Each ligand was docked in three runs with seeds 1, 25 and 78 using AutoDock Vina 1.2.5.

## Molecular dynamics

System labels: apo FADS1; C1 = FADS1–Luminespib; C2 = FADS1–NU7441; T = FADS1–D5D-IN-326. Each system was simulated once for 100 ns in water (CHARMM36, TIP3P, 298 K, 1 bar) with GROMACS 2026.3.

- `mdp_files/`: `ions.mdp`, `minim.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp` (10 ns), `md_50ns.mdp`, `md_100ns.mdp`, `ie.mdp` (interaction-energy rerun) and `mmpbsa.in` (MM-GBSA, frames 8,001–10,000 every 10 frames, igb = 5, 0.150 M salt). The run scripts expect `mdp_files/` under `${MD_ROOT}/fads1_md_c1/` (C1, C2) and `${MD_ROOT}/fads1_md/` (apo, T), and `ie.mdp` and `mmpbsa.in` in each system folder.
- `setup/`: apo and T system build, equilibration and production (`setup_apo_and_d5d_in_326.sh`), extension of the apo run from 10 to 100 ns (`extend_apo_to_100ns.sbatch`), C1 equilibration and production (`run_luminespib.sbatch`), C2 complex build, equilibration, 50-ns production and extension to 100 ns (`nu7441_*.sbatch`, `merge_gro.py`, `check_complex.py`).
- `analysis/`: trajectory analysis (RMSD, RMSF, radius of gyration, hydrogen bonds, SASA, interaction energy, PCA) for C1 (`analysis_luminespib.sh`), T (`analysis_d5d_in_326.sh`) and apo and C2 (`postmd_analysis.sbatch`), and MM-GBSA (`run_mmgbsa.sh` and wrappers).

SLURM directives reflect the University of Tartu HPC cluster and should be adapted to other systems.

## Software

See `environment/`. Main versions: R 4.5.2 with TwoSampleMR 0.7.4, MVMR 0.4.4, MRPRESSO 1.0, MendelianRandomization 0.10.0, mr.raps 0.4.3, ieugwasr 1.1.0, coloc 5.2.3, susieR 0.14.2 and data.table 1.18.4; LDSC; S-PrediXcan; PLINK 1.9; AutoDock Vina 1.2.5; Meeko 0.7.1; GROMACS 2026.3; ACPYPE 2026.9.4; gmx_MMPBSA 1.7.0.

## License

MIT (see `LICENSE`). If you use this code, please cite the article (see `CITATION.cff`).
