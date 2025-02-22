##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


library(enrichR)
require(httr)

listEnrichrSites()

dbs <- listEnrichrDbs()

# remove kegg
dbs <- dbs[grep("KEGG", dbs$libraryName, invert = TRUE), ]
dbs <- dbs[grep("kegg", dbs$link, invert = TRUE), ]

# remove mouse (we are only interested in human pathways, mouse should be added in the future as a separate species)
dbs <- dbs[grep("mouse", dbs$link, invert = TRUE), ]

# remove wikipathways, as enrichr excludes some pathways and we dont want that (we want full WP)
dbs <- dbs[grep("WikiPathways", dbs$libraryName, invert = TRUE), ]
dbs <- dbs[grep("WikiPathway", dbs$libraryName, invert = TRUE), ]

# manually remove cases where the numeric year is not the last element

dbs <- dbs[!dbs$libraryName %in% c(
    "GO_Biological_Process_2017b",
    "GO_Cellular_Component_2017b",
    "GO_Molecular_Function_2017b"
), ]

# Get the year when possible

split_string <- strsplit(dbs$libraryName, "_")

last_element <- lapply(split_string, function(x) {
    year <- tail(x, n = 1)
    year <- ifelse(grepl("^[0-9]+$", year), as.numeric(year), NA)
})

last_element <- unlist(last_element)

dbs$year <- last_element

# remove the year from the library name
dbs$DS_Name <- gsub("_[0-9]+$", "", dbs$libraryName)
# keep the latest library according to the year, if there's more than one library with the same name
dbs <- dbs[order(dbs$year, decreasing = TRUE), ]
dbs <- dbs[!duplicated(dbs$DS_Name), ]


# api call to retrieve the gene sets

lapply(dbs$libraryName, function(db_name) {
    # db_name <- dbs$libraryName[1]
    url <- paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", db_name)
    resp <- GET(url)
    filename <- paste0(db_name, "_AUTOGENERATED", ".txt")

    write(rawToChar(resp$content), file = file.path("data-raw", "extdata", "gmt", filename))
})

write.csv(dbs, file.path("data-raw", "extdata", "gmt", "info_enrichR.csv"), row.names = FALSE)
