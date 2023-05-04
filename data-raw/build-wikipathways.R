##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

if(!"rWikiPathways" %in% installed.packages()){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("rWikiPathways", update = FALSE)
}

require(rWikiPathways)
require(biomaRt)

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Read wikipathways GMT file


# eventually we can download wikipathways from different species with 
# https://bioconductor.org/packages/release/bioc/html/rWikiPathways.html
gmt <- readGMT('data-raw/extdata/wikipathways-20230410-gmt-Homo_sapiens.gmt')

# Retrieve the gene symbols using the getBM function

gene_symbols = getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                     filters = "entrezgene_id",
                     values = gmt[2],
                     mart = ensembl)
gmt$gene <- gene_symbols[match(gmt[,2], gene_symbols[,1]),2]

writeGMT(gmt, "data-raw/extdata/gmt/Pathway_wikipathways_hsa.gmt")

