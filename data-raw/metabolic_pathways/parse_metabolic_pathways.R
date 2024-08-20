library(dplyr)

# Read file
data <- readr::read_tsv("./data-raw/metabolic_pathways/ChEBI2Reactome.txt")

# Inspect the data to determine the correct number of columns
head(data)

colnames(data) <- c("chebi", "stId", "reactome_url", "reactome_name", "evidence", "species")

head(data)

# keep Homo sapiens for now
data <- data[data$species == "Homo sapiens", ]

head(data)
dim(data)

# append reactome path is to pathways name
data$reactome_name <- paste(data$reactome_name, data$stId, sep = "_")


gmt_data <- data %>%
    group_by(reactome_name) %>%
    summarise(chebi_ids = paste(chebi, collapse = "\t")) %>%
    ungroup()


write.table(gmt_data,
    file = "./data-raw/metabolic_pathways/metabolic_pathways.gmt",
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# clean up genesets
gmt <- playbase::read.gmt(gmt.file = "./data-raw/metabolic_pathways/metabolic_pathways.gmt")
gmt <- list(CUSTOM = gmt)
gmt <- playbase::clean_gmt(gmt, "PATHWAY")

# compute custom geneset stats
gmt <- gmt[!duplicated(names(gmt))]
gset_size <- sapply(gmt, length)
gmt <- gmt[gset_size >= 3]

REACTOME_METABOLITES <- gmt
usethis::use_data(REACTOME_METABOLITES, overwrite = TRUE)


# parse Wikipathwaays
library(AnnotationHub)

ah <- AnnotationHub()

query(ah, "wikipathways")

qr <- query(ah, c("wikipathways", "Homo sapiens"))

hsatbl <- qr[[1]]

colnames(hsatbl)

head(hsatbl)

sum(duplicated(hsatbl$ChEBI_ID))

length(unique(hsatbl$ChEBI_ID))

length(unique(hsatbl$wpid))

data <- hsatbl

# append reactome path is to pathways name

# keep only numerical characters from ChEBI_ID
data$ChEBI_ID <- gsub("[^0-9]", "", data$ChEBI_ID)

data$wp_name <- paste(data$pathway_name, data$wpid, sep = "_")

gmt_data <- data %>%
    group_by(wp_name) %>%
    summarise(chebi_ids = paste(ChEBI_ID, collapse = "\t")) %>%
    ungroup()


write.table(gmt_data,
    file = "./data-raw/metabolic_pathways/metabolic_pathways_wp.gmt",
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# clean up genesets
gmt <- playbase::read.gmt(gmt.file = "./data-raw/metabolic_pathways/metabolic_pathways_wp.gmt")
gmt <- list(CUSTOM = gmt)
gmt <- playbase::clean_gmt(gmt, "PATHWAY")

# compute custom geneset stats
gmt <- gmt[!duplicated(names(gmt))]
gset_size <- sapply(gmt, length)
gmt <- gmt[gset_size >= 3]

WP_METABOLITES <- gmt
usethis::use_data(WP_METABOLITES, overwrite = TRUE)
