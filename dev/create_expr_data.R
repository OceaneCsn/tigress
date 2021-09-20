library(DIANE)
data("abiotic_stresses")
data("gene_annotations")
data("regulators_per_organism")

# getting raw counts from DIANE's demo data
tcc_object <- list(counts = abiotic_stresses$raw_counts)
threshold = 10*length(abiotic_stresses$conditions)
tcc_object$counts <- tcc_object$counts[rowSums(tcc_object$counts) > threshold,]
normalized_counts <- tcc_object$counts

fit <- DIANE::estimateDispersion(tcc = tcc_object)


topTags <- DIANE::estimateDEGs(fit, reference = "M",
                               perturbation = "MH",
                               p.value = 0.05, lfc = 1.5)

genes <- get_locus(topTags$table$genes)
regressors <- intersect(genes, regulators_per_organism[[
  "Arabidopsis thaliana"]])

aggregated_data <- aggregate_splice_variants(data = normalized_counts)

grouping <- DIANE::group_regressors(aggregated_data, genes,
                                    regressors, corr_thr = 0.95)

grouped_counts <- grouping$counts
grouped_targets <- grouping$grouped_genes
grouped_regressors <- grouping$grouped_regressors


abiotic_stresses_data <- list(expression = grouped_counts,
                              genes = grouped_targets,
                              regulators = grouped_regressors)
