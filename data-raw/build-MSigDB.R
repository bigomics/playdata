##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


library(msigdbr)

all_species <- msigdbr_species()

all_gene_sets_xt = msigdbr(species = all_species[20,1][[1]])

all_gene_sets = msigdbr(species = "Homo sapiens")

gmt <-  playbase::convert.gmt(
    all_gene_sets$gene_symbol,
    gs_name = all_gene_sets$gs_name)

# add url as first gene in case its necessary
# gmt_2 <- lapply(names(gmt), function(x){
#     genes <- gmt[[x]]
#     url <- unique(url[[x]])
#     if(url ==""){
#         return(genes)
#     }
#     return(c(url, genes))
# })
# url <- playbase::convert.gmt(
#     all_gene_sets$gs_url,
#     gs_name = all_gene_sets$gs_name)
# names(gmt_2) <- names(gmt)
# gmt_2["ADA2_TARGET_GENES"]

playbase::write.gmt(gmt, file= file.path("data-raw","extdata","gmt_msigdb","all_gene_sets.txt"))
