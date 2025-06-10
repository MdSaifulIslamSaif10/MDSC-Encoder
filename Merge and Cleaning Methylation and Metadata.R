###### Load DNA Methylation Data
load("R:/BIO RSCH/BioStats/EPIC_data/Sudhi/Betas.RData")
dim(Betas)
colnames(Betas)
rownames(Betas)

Betas_t <- t(Betas)

###### Load Metadata
metadata <- read.csv("R:/BIO RSCH/BioStats/EPIC_data/Metadata for Saif from Dr Koestler/metadata.csv")
SampleId_for_merge <- read.csv("R:/BIO RSCH/BioStats/EPIC_data/Metadata for Saif from Dr Koestler/sample_manifest.csv")

metadata$Sample_ID
SampleId_for_merge$Sample_Name

###### Merge Metadata
merged_metadata_modified_id <- merge(SampleId_for_merge, metadata,
                                     by.x = "Sample_Name", by.y = "Sample_ID",
                                     all.x = TRUE)

sum(is.na(merged_metadata_modified_id$X..M.MDSCs_PBMCs))
sum(is.na(merged_metadata_modified_id$X..G.MDSCs_PBMCs))

###### Prepare for Merging with Methylation Data
merged_metadata_modified_id$methylation_id_for_merge <- paste0(
  "X",
  merged_metadata_modified_id$Sentrix_ID,
  "_",
  merged_metadata_modified_id$Sentrix_Position
)

Betas_t <- as.data.frame(Betas_t)
Betas_t <- cbind(methylation_id_for_merge = rownames(Betas_t), Betas_t)
rownames(Betas_t) <- NULL

Methylation_and_metadata <- merge(Betas_t, merged_metadata_modified_id,
                                  by = "methylation_id_for_merge", all.x = TRUE)

###### Extract and Clean Beta Matrix
beta_raw <- as.matrix(Methylation_and_metadata[, 2:936991])
rownames(beta_raw) <- Methylation_and_metadata$methylation_id_for_merge
clean_cpg_ids <- sub("_.*", "", colnames(beta_raw))
colnames(beta_raw) <- clean_cpg_ids
beta_matrix <- t(beta_raw)

###### Load EPIC Annotation and Filter Probes
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

sex_probes <- rownames(annotation[annotation$chr %in% c("chrX", "chrY"), ])
snp_probes <- rownames(annotation[!is.na(annotation$Probe_rs), ])
probes_to_remove <- unique(c(sex_probes, snp_probes))

beta_filtered <- beta_matrix[!(rownames(beta_matrix) %in% probes_to_remove), ]

dim(beta_filtered)
head(rownames(beta_filtered))

###### Remove Cross-Reactive Probes
library(maxprobes)
cross_reactive_probes <- maxprobes::xreactive_probes(array_type = "EPIC")
beta_final <- beta_filtered[!(rownames(beta_filtered) %in% cross_reactive_probes), ]

dim(beta_final)
head(rownames(beta_final))

###### Filter Merged Data to Keep Only Good Probes
methylation_data <- Methylation_and_metadata[, 2:936991]
cleaned_names <- sub("_.*", "", colnames(methylation_data))
colnames(methylation_data) <- cleaned_names
good_idx <- which(cleaned_names %in% rownames(beta_final))
filtered_methylation_data <- methylation_data[, good_idx]

meta_idx <- setdiff(seq_len(ncol(Methylation_and_metadata)), 2:936991)
metadata <- Methylation_and_metadata[, meta_idx]

Methylation_and_metadata_filtered <- cbind(metadata, filtered_methylation_data)

dim(Methylation_and_metadata_filtered)
head(colnames(Methylation_and_metadata_filtered))

###### Final Subsetting to CpGs Only
cpg_cols <- colnames(Methylation_and_metadata_filtered)[55:799688]
cg_only_cols <- cpg_cols[grepl("^cg", cpg_cols)]
metadata_cols <- colnames(Methylation_and_metadata_filtered)[1:54]

Methylation_and_metadata_filtered <- Methylation_and_metadata_filtered[, c(metadata_cols, cg_only_cols)]

dim(Methylation_and_metadata_filtered)
head(colnames(Methylation_and_metadata_filtered))

###### Select Top Variable CpGs
metadata_cols <- 1:54
cpg_cols <- (max(metadata_cols) + 1):ncol(Methylation_and_metadata_filtered)

cpg_sd <- apply(Methylation_and_metadata_filtered[, cpg_cols], 2, sd, na.rm = TRUE)
n_top <- ceiling(.01 * length(cpg_sd))
top_cpg_indices_within <- order(cpg_sd, decreasing = TRUE)[1:n_top]
top_cpg_cols <- cpg_cols[top_cpg_indices_within]

filtered_data <- Methylation_and_metadata_filtered[, c(metadata_cols, top_cpg_cols)]
dim(filtered_data)
