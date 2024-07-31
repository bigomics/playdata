library(dplyr)

# Read file
data <- read.delim("./data-raw/metabolic_pathways/ChEBI2Reactome.txt", header = FALSE)

# Inspect the data to determine the correct number of columns
head(data)

colnames(data) <- c("chebi", "stId", "reactome_url", "reactome_name", "evidence", "species")

head(data)

gmt_data <- data %>%
    group_by(reactome_name, reactome_url) %>%
    summarise(chebi_ids = paste(chebi, collapse = "\t")) %>%
    ungroup()


write.table(gmt_data,
    file = "./data-raw/metabolic_pathways.gmt",
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)
