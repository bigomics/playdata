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

playbase::write.gmt(gmt, file= file.path("data-raw","extdata","gmt_msigdb","all_gene_sets.txt"))

meta <- all_gene_sets[!duplicated(all_gene_sets[ , c("gs_name")]),]
write.csv(meta ,file.path("data-raw","extdata","gmt_msigdb","meta.csv"))