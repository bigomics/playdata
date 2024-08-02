library(dplyr)

# Read file
data_chebi <- read.delim("./data-raw/metabolic_pathways/chebi_compounds_20240801_0501.tsv", header = TRUE)

head(data_chebi)

# remove unnecessary columns
data_chebi <- data_chebi[c("ID", "CHEBI_ACCESSION", "NAME", "DEFINITION")]

data_conversion <- read.csv("./data-raw/metabolic_pathways/mapping_metaboanalyst_20240802_1508.csv", header = TRUE)

dim(data_conversion)
# remove comment column by name comment
data_conversion <- data_conversion[, !(names(data_conversion) %in% c("Comment", "SMILES"))]

# convert all columns to character
data_conversion[] <- lapply(data_conversion, as.character)

head(data_conversion)

dim(data_conversion)

# remove rows that are all na, except for the first column
data_conversion <- data_conversion[rowSums(is.na(data_conversion[, -1])) < ncol(data_conversion) - 1, ]

# count unique ids
length(unique(data_conversion$Query))

# merge the two dataframes
METABOLITE_ANNOTATION <- merge(data_chebi, data_conversion, by.x = "ID", by.y = "Query", all.x = TRUE)

dim(METABOLITE_ANNOTATION)

unique(METABOLITE_ANNOTATION$KEGG)

colnames(METABOLITE_ANNOTATION)

# check how many ids map to metabolic pathways
reactome <- playdata::REACTOME_METABOLITES

# get the unique ids from all lists
reactome_ids <- unique(unlist(reactome))
length(reactome_ids)

sum(METABOLITE_ANNOTATION$ID %in% reactome_ids)

# subset

metab_in_pathway <- METABOLITE_ANNOTATION[METABOLITE_ANNOTATION$ID %in% reactome_ids, ]

dim(metab_in_pathway)

unique(metab_in_pathway$KEGG)
unique(metab_in_pathway$HMDB)

usethis::use_data(METABOLITE_ANNOTATION, overwrite = TRUE)
