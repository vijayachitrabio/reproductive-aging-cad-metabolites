#!/usr/bin/env Rscript
# STRING physical-interaction network for the FADS1 seed set: centrality and
# Louvain modules (Figure 5, Supplementary Tables S11-S12).

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

root <- normalizePath(Sys.getenv("RAC_ROOT", "."), mustWork = FALSE)
if (!dir.exists(file.path(root, "results", "ppi_network"))) dir.create(file.path(root, "results", "ppi_network"), recursive = TRUE)
out <- file.path(root, "results", "ppi_network")
fig <- file.path(root, "manuscript", "figures")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

seed <- read_tsv(file.path(out, "seed_genes.tsv"), show_col_types = FALSE)
net <- read_tsv(file.path(out, "string_network_expanded.tsv"), show_col_types = FALSE)
net <- net %>% transmute(from = preferredName_A, to = preferredName_B, score = as.numeric(score)) %>%
  filter(!is.na(from), !is.na(to), from != to) %>% distinct()

g <- graph_from_data_frame(net, directed = FALSE)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(score = "max"))
V(g)$seed <- V(g)$name %in% seed$gene
V(g)$tier <- seed$evidence_tier[match(V(g)$name, seed$gene)]
V(g)$tier[is.na(V(g)$tier)] <- "added_STRING_interactor"
E(g)$weight <- E(g)$score

harm <- if ("harmonic_centrality" %in% getNamespaceExports("igraph")) harmonic_centrality(g, weights = NULL) else closeness(g, normalized = TRUE)
metrics <- tibble(
  gene = V(g)$name,
  seed = V(g)$seed,
  evidence_tier = V(g)$tier,
  degree = degree(g),
  betweenness = betweenness(g, directed = FALSE, normalized = TRUE),
  closeness_or_harmonic = as.numeric(harm),
  eigenvector = eigen_centrality(g, weights = E(g)$weight, directed = FALSE)$vector,
  kcore = coreness(g)
)
metric_cols <- c("degree", "betweenness", "closeness_or_harmonic", "eigenvector", "kcore")
for (m in metric_cols) metrics[[paste0(m, "_percentile")]] <- percent_rank(metrics[[m]])
metrics <- metrics %>%
  mutate(consensus_score = rowMeans(across(ends_with("_percentile")), na.rm = TRUE),
         top_consensus = consensus_score >= quantile(consensus_score, 0.80, na.rm = TRUE),
         bottleneck = betweenness_percentile >= quantile(betweenness_percentile, 0.80, na.rm = TRUE),
         hub = rowSums(across(all_of(paste0(metric_cols, "_percentile"))) >= 0.80) >= 3,
         hub_and_bottleneck = hub & bottleneck) %>% arrange(desc(consensus_score))
write_tsv(metrics, file.path(out, "network_centrality_metrics.tsv"))

set.seed(260926)
comm <- cluster_louvain(g, weights = E(g)$weight)
membership_tbl <- tibble(gene = names(membership(comm)), module = as.integer(membership(comm))) %>%
  left_join(metrics %>% select(gene, seed, evidence_tier), by = "gene")
write_tsv(membership_tbl, file.path(out, "network_module_membership.tsv"))
write_tsv(net, file.path(out, "string_network_edges.tsv"))

fads_nodes <- ego(g, order = 1, nodes = "FADS1", mindist = 0)[[1]]
fads_edges <- igraph::as_data_frame(induced_subgraph(g, vids = fads_nodes), what = "edges")
write_tsv(fads_edges, file.path(out, "fads1_first_neighborhood_edges.tsv"))

theme_set(theme_minimal(base_size = 11))
top <- metrics %>% filter(seed) %>% slice_max(consensus_score, n = min(15, nrow(.))) %>%
  mutate(gene = reorder(gene, consensus_score))
p1 <- ggplot(top, aes(consensus_score, gene, colour = hub_and_bottleneck, shape = evidence_tier)) +
  geom_point(size = 3) + labs(title = "Consensus centrality among seeded genes", x = "Mean centrality percentile", y = NULL, colour = "Hub + bottleneck") +
  scale_colour_manual(values = c(`FALSE` = "#4C78A8", `TRUE` = "#D1495B"))
ggsave(file.path(fig, "Figure_PPI_consensus_centrality.png"), p1, width = 8, height = 5.4, dpi = 600)

lay <- layout_with_fr(g, weights = E(g)$weight, niter = 2000)
png(file.path(fig, "Figure_PPI_network_overview.png"), width = 2600, height = 2200, res = 350)
plot(g, layout = lay, vertex.label = ifelse(V(g)$seed, V(g)$name, NA), vertex.size = ifelse(V(g)$name == "FADS1", 13, ifelse(V(g)$seed, 8, 5)), vertex.color = ifelse(V(g)$name == "FADS1", "#D1495B", ifelse(V(g)$seed, "#4C78A8", "#BDBDBD")), vertex.frame.color = "white", vertex.label.cex = 0.65, vertex.label.color = "#222222", edge.width = 1 + 3 * E(g)$weight, edge.color = "#B8B8B8", main = "STRING physical interaction network: ovarian-aging/FADS1 exploratory seed set")
legend("topleft", legend = c("FADS1", "seed gene", "added STRING interactor"), pch = 21, pt.cex = c(1.7, 1.3, 0.9), pt.bg = c("#D1495B", "#4C78A8", "#BDBDBD"), bty = "n")
dev.off()

sg <- induced_subgraph(g, vids = fads_nodes)
png(file.path(fig, "Figure_FADS1_centered_subnetwork.png"), width = 2600, height = 2200, res = 350)
plot(sg, layout = layout_with_fr(sg, weights = E(sg)$weight, niter = 1500), vertex.label = V(sg)$name, vertex.size = ifelse(V(sg)$name == "FADS1", 16, 10), vertex.color = ifelse(V(sg)$name == "FADS1", "#D1495B", ifelse(V(sg)$seed, "#4C78A8", "#BDBDBD")), vertex.frame.color = "white", vertex.label.cex = 0.8, edge.width = 1 + 3 * E(sg)$weight, edge.color = "#A9A9A9", main = "FADS1 first-neighborhood physical interaction network")
dev.off()

writeLines(c(
  paste0("nodes\t", vcount(g)), paste0("edges\t", ecount(g)), paste0("seed_nodes\t", sum(V(g)$seed)),
  paste0("added_STRING_interactors\t", sum(!V(g)$seed)), paste0("modules\t", length(unique(membership(comm)))),
  paste0("FADS1_degree\t", metrics$degree[metrics$gene == "FADS1"]),
  paste0("FADS1_betweenness\t", signif(metrics$betweenness[metrics$gene == "FADS1"], 5)),
  paste0("FADS1_consensus_percentile\t", signif(metrics$consensus_score[metrics$gene == "FADS1"], 5)),
  paste0("FADS1_hub\t", metrics$hub[metrics$gene == "FADS1"]),
  paste0("FADS1_bottleneck\t", metrics$bottleneck[metrics$gene == "FADS1"])
), file.path(out, "network_summary.tsv"))

message("Wrote network metrics, module membership, and figures to ", out)
