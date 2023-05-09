##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


library(msigdbr)

all_species <- msigdbr_species()

all_gene_sets_xt = msigdbr(species = all_species[20,1][[1]])

all_gene_sets = msigdbr(species = "Homo sapiens")

dim(all_gene_sets_xt)

dim(all_gene_sets)
