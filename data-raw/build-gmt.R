##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


# script to create GSETS used in omicsplayground

# get wikipathways genesets
source("data-raw//build-gmt//get-wikipathways.R")

# get single databases genesets
source("data-raw//build-gmt//get-gmt.r")

# get MSigDB genesets
source("data-raw//build-gmt//get-MSigDB.R")

# build gmt-all file, that merges all genesets above
source("data-raw//build-gmt//build-gmt-all.R")

# create sparse XL matrix with genesets as rows and genes as columns
source("data-raw//build-gmt//build-GSET_SPARSEG_XL.R")

# build families
source("data-raw//build-data.R")

