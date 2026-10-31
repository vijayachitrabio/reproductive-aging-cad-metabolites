#!/usr/bin/env Rscript
# Network robustness across STRING settings and evidence overlay for network genes
# (Table 4, Supplementary Tables S10-S13, S15).

suppressPackageStartupMessages({
  library(igraph)
  library(readr)
  library(dplyr)
  library(tidyr)
})

root <- normalizePath(Sys.getenv("RAC_ROOT", "."), mustWork = FALSE)
netdir <- file.path(root, "results", "ppi_network")
robdir <- file.path(netdir, "robustness")
evidir <- file.path(netdir, "evidence_overlay")
dir.create(robdir, recursive = TRUE, showWarnings = FALSE)
dir.create(evidir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  baseline_score700_add10 = file.path(robdir, "string_network_score700_add10_baseline.tsv"),
  high_score800_add5 = file.path(robdir, "string_network_score800_add5.tsv"),
  high_score800_add0 = file.path(robdir, "string_network_score800_add0.tsv")
)

read_net <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(from = preferredName_A, to = preferredName_B, score = as.numeric(score)) %>%
    filter(!is.na(from), !is.na(to), from != to) %>% distinct()
}

seed <- read_tsv(file.path(netdir, "seed_genes.tsv"), show_col_types = FALSE)

analyze_one <- function(label, path) {
  e <- read_net(path)
  g <- simplify(graph_from_data_frame(e, directed = FALSE), remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(score = "max"))
  E(g)$weight <- E(g)$score
  V(g)$seed <- V(g)$name %in% seed$gene
  if (vcount(g) > 2 && ecount(g) > 1) comm <- cluster_louvain(g, weights = E(g)$weight) else comm <- make_clusters(g, membership = seq_len(vcount(g)))
  mem <- membership(comm)
  fads_module <- if ("FADS1" %in% V(g)$name) unname(mem["FADS1"]) else NA_integer_
  fads_members <- if (!is.na(fads_module)) names(mem[mem == fads_module]) else character()
  harm <- if ("harmonic_centrality" %in% getNamespaceExports("igraph")) harmonic_centrality(g, weights = NULL) else closeness(g, normalized = TRUE)
  metrics <- tibble(
    setting = label,
    gene = V(g)$name,
    seed = V(g)$seed,
    degree = degree(g),
    betweenness = betweenness(g, directed = FALSE, normalized = TRUE),
    closeness_or_harmonic = as.numeric(harm),
    eigenvector = eigen_centrality(g, weights = E(g)$weight, directed = FALSE)$vector,
    module = as.integer(mem[V(g)$name])
  )
  fads <- metrics %>% filter(gene == "FADS1")
  summary <- tibble(
    setting = label,
    threshold = ifelse(grepl("700", label), 700, 800),
    added_interactors = ifelse(grepl("add10", label), 10, ifelse(grepl("add5", label), 5, 0)),
    nodes = vcount(g), edges = ecount(g), modules = length(unique(mem)),
    observed_seed_nodes = sum(V(g)$seed), absent_input_seeds = sum(!seed$gene %in% V(g)$name),
    fads1_degree = ifelse(nrow(fads), fads$degree, NA),
    fads1_betweenness = ifelse(nrow(fads), fads$betweenness, NA),
    fads1_module = fads_module,
    fads1_module_size = length(fads_members),
    fads1_module_members = paste(sort(fads_members), collapse = ",")
  )
  modules <- tibble(setting = label, gene = names(mem), module = as.integer(mem)) %>%
    group_by(setting, module) %>% mutate(module_size = n()) %>% ungroup()
  list(summary = summary, metrics = metrics, modules = modules, graph = g)
}

results <- Map(analyze_one, names(files), unname(files))
summaries <- bind_rows(lapply(results, `[[`, "summary"))
metrics <- bind_rows(lapply(results, `[[`, "metrics"))
modules <- bind_rows(lapply(results, `[[`, "modules"))
write_tsv(summaries, file.path(robdir, "robustness_summary.tsv"))
write_tsv(metrics, file.path(robdir, "robustness_centrality_metrics.tsv"))
write_tsv(modules, file.path(robdir, "robustness_module_membership.tsv"))

module_sets <- lapply(results, function(x) {
  s <- x$summary$fads1_module_members
  if (!nzchar(s)) character() else strsplit(s, ",", fixed = TRUE)[[1]]
})
jaccard <- function(a, b) if (!length(a) && !length(b)) 1 else if (!length(a) || !length(b)) 0 else length(intersect(a, b)) / length(union(a, b))
agreement <- expand.grid(setting_a = names(files), setting_b = names(files), stringsAsFactors = FALSE) %>%
  rowwise() %>% mutate(jaccard_fads1_module = jaccard(module_sets[[setting_a]], module_sets[[setting_b]])) %>% ungroup()
write_tsv(agreement, file.path(robdir, "fads1_module_agreement_jaccard.tsv"))

# Evidence overlay: direct project evidence is kept separate from curated pathway context.
all_genes <- sort(unique(unlist(lapply(results, function(x) V(x$graph)$name))))
fatty_context <- c("FADS2", "ELOVL2", "ELOVL5", "SCD", "ACOX1", "PLA2G4A", "PTGS1", "PTGS2", "ALOX5", "ALOX15")
ovarian_context <- c("AMH", "AMHR2", "FOXL2", "FST", "FSHR", "LHCGR", "STAR", "CYP19A1", "CYP17A1", "HSD17B1", "ESR1")
overlay <- tibble(
  gene = all_genes,
  anm_genetic_evidence = gene %in% "FADS1",
  amh_genetic_evidence = FALSE,
  fads1_cis_eqtl_or_coloc = gene %in% "FADS1",
  fatty_acid_or_metabolite_evidence = gene %in% "FADS1",
  ovarian_bulk_transcriptomic_evidence = gene %in% c("FADS1", "FADS2"),
  granulosa_specific_transcriptomic_evidence = FALSE,
  curated_fatty_acid_pathway_context = gene %in% fatty_context,
  curated_ovarian_follicular_context = gene %in% ovarian_context,
  string_only_added_interactor = !(gene %in% seed$gene)
) %>%
  mutate(independent_evidence_types = rowSums(across(c(anm_genetic_evidence, amh_genetic_evidence, fads1_cis_eqtl_or_coloc, fatty_acid_or_metabolite_evidence, ovarian_bulk_transcriptomic_evidence))),
         multiple_independent_sources = independent_evidence_types >= 2,
         evidence_note = case_when(
           gene == "FADS1" ~ "ANM-Omega3 coloc PP4 approximately 0.71; FADS1 cis-eQTL/coloc context; fatty-acid metabolite prioritization; bulk ovary RNA support; not granulosa-specific",
           gene == "FADS2" ~ "bulk ovary RNA support plus curated FADS pathway context; no gene-level ANM/AMH evidence identified in project tables",
           string_only_added_interactor ~ "Added by STRING high-confidence physical-network expansion; no independent project evidence assigned",
           gene %in% fatty_context ~ "Curated fatty-acid/eicosanoid pathway context; not a direct project gene-level association",
           gene %in% ovarian_context ~ "Curated ovarian/follicular context; not a direct gene-level AMH or granulosa result in available project tables",
           TRUE ~ "No independent project evidence identified"
         ))
write_tsv(overlay, file.path(evidir, "network_gene_evidence_overlay.tsv"))

module_evidence <- modules %>% left_join(overlay, by = "gene") %>%
  group_by(setting, module, module_size) %>%
  summarise(genes = paste(sort(gene), collapse = ","),
            independent_source_types_present = paste(c("ANM_genetic", "AMH_genetic", "FADS1_cis_eqtl_or_coloc", "fatty_acid_or_metabolite", "ovarian_bulk_transcriptomic")[c(any(anm_genetic_evidence), any(amh_genetic_evidence), any(fads1_cis_eqtl_or_coloc), any(fatty_acid_or_metabolite_evidence), any(ovarian_bulk_transcriptomic_evidence))], collapse = ";"),
            genes_with_multiple_independent_sources = paste(sort(gene[multiple_independent_sources]), collapse = ","),
            .groups = "drop") %>%
  mutate(n_independent_source_types = ifelse(independent_source_types_present == "", 0, lengths(strsplit(independent_source_types_present, ";", fixed = TRUE))),
         contains_fads1 = grepl("(^|,)FADS1(,|$)", genes))
write_tsv(module_evidence, file.path(evidir, "module_evidence_summary.tsv"))

classify <- function(x, fads_members) case_when(
  x %in% fads_members ~ "FADS1_lipid_module",
  x %in% c("AMH", "AMHR2") ~ "AMH_AMHR2_module",
  x %in% c("FSHR", "LHCGR", "FSHB", "CGA") ~ "follicular_module",
  x %in% c("PTGS1", "PTGS2") ~ "cyclooxygenase_module",
  x %in% c("ALOX5", "ALOX5AP", "ALOX15", "ALOX15B") ~ "lipoxygenase_module",
  TRUE ~ "other"
)
inter_edges <- bind_rows(lapply(results, function(x) {
  fads_members <- strsplit(x$summary$fads1_module_members, ",", fixed = TRUE)[[1]]
  if (!length(fads_members) || !nzchar(fads_members[1])) fads_members <- character()
  igraph::as_data_frame(x$graph, what = "edges") %>% mutate(
    setting = x$summary$setting,
    class_from = classify(from, fads_members),
    class_to = classify(to, fads_members)
  ) %>% filter(class_from != class_to, class_from != "other", class_to != "other")
}))
write_tsv(inter_edges, file.path(robdir, "inter_module_edges_all_settings.tsv"))
write_tsv(inter_edges %>% filter(setting == "baseline_score700_add10"), file.path(robdir, "baseline_inter_module_edges.tsv"))

message("Wrote robustness, overlay and inter-module outputs to ", netdir)
