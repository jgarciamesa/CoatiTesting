# Get gene ID list from Ensembl

library(biomaRt)
library(dplyr, warn.conflicts = FALSE)

# As explained in bioconductor.org, if R version is 3.6 or above:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("biomaRt")
#
# For documentation of this version:
# browseVignettes("biomaRT")

dir.create(file.path("raw_data","dmelanogaster"), showWarnings = FALSE)

#output = "raw_data/mouse/mouse_geneId.tsv"
output = "raw_data/dmelanogaster_geneId.tsv"

# Query to Ensembl and select homo sapiens gene dataset
ensembl_h = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Filters to use
filters = c("with_ccds")
# Attributes to retrieve
#attributes = c("ensembl_gene_id", "ggorilla_homolog_ensembl_gene", "ggorilla_homolog_orthology_type")
#attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type")
attributes = c("ensembl_gene_id", "dmelanogaster_homolog_ensembl_gene", "dmelanogaster_homolog_orthology_type")
values = list(TRUE)#, TRUE)

human_otherspecies = getBM(attributes, filters, values,
	mart = ensembl_h, uniqueRows = TRUE)

# Filter one one2one (121) orthologs
human_otherspecies_121_orth = filter(human_otherspecies, dmelanogaster_homolog_orthology_type == "ortholog_one2one")

write.table(human_otherspecies_121_orth[,1:2], file = output, quote = FALSE,
	sep = "\t", row.names = FALSE, col.names = TRUE)

