## NOTE: see original file global-init.R (not to be used anymore)
## NOTE2: please check if these are still actively being used in OPG or playbase
## NOTE3: we have replaced 'dots' with 'underscore'! please note.

library(playbase)

## These source files are still in OPG but probably be moved here.
FILES = "~/Playground/omicsplayground/lib"

##---------------------------------------------------------
## All gene families in Human UPPER CASE
require(org.Hs.eg.db)
GENE_TITLE  = unlist(as.list(org.Hs.egGENENAME))
GENE_SYMBOL = unlist(as.list(org.Hs.egSYMBOL))
names(GENE_TITLE) = GENE_SYMBOL
usethis::use_data(GENE_TITLE, overwrite = TRUE)
usethis::use_data(GENE_SYMBOL, overwrite = TRUE)

##GSET.PREFIX.REGEX = paste(paste0("^",GSET.PREFIXES,"_"),collapse="|")
GSET_PREFIX_REGEX="^BIOCARTA_|^C2_|^C3_|^C7_|^CHEA_|^GOBP_|^GOCC_|^GOMF_|^HALLMARK_|^KEA_|^KEGG_|^PID_|^REACTOME_|^ST_"
usethis::use_data(GSET_PREFIX_REGEX)

##---------------------------------------------------------
GENE_SUMMARY = read.csv(file.path(FILES,"gene-summary.csv"),row.names=1)
GENE_SUMMARY = array(GENE_SUMMARY[,1], dimnames=list(rownames(GENE_SUMMARY)))
usethis::use_data(GENE_SUMMARY)

##---------------------------------------------------------
## GENExGENE <- readRDS(file=file.path(FILES,"GENExGENE-cosSparseKNN500-XL.rds"))
GSETxGENE <- readRDS(file.path(FILES,"gset-sparseG-XL.rds"))
usethis::use_data(GSETxGENE)

##---------------------------------------------------------
load(file.path(FILES,"gmt-all.rda"),verbose=1)
GSETS = gmt.all;remove(gmt.all)
##usethis::use_data(GSETS)

##---------------------------------------------------------
message("[INIT] parsing gene families...")
FAMILIES <- pgx.getGeneFamilies(GENE_SYMBOL, FILES=FILES, min.size=10, max.size=9999)
#fam.file <- file.path(FILES,"custom-families.gmt")
#if(file.exists(fam.file)) {
#  custom.gmt = read.gmt(file.path(FILES,"custom-families.gmt"),add.source=TRUE)
#  names(custom.gmt)
#  FAMILIES= c(FAMILIES, custom.gmt)
#}
FAMILIES[["<all>"]] <- GENE_SYMBOL
f1 <- FAMILIES
names(f1) <- paste0("FAMILY:",names(f1))
names(f1) <- sub("FAMILY:<all>","<all>",names(f1))
GSETS <- c(GSETS,f1)

usethis::use_data(FAMILIES)
usethis::use_data(GSETS, overwrite=TRUE)

##---------------------------------------------------------
## The GSETS object is really really big because its a big list of
## gene names. It is more efficient to store indices/integers to a
## gene list. iGSETS is a list of integers and needs much smaller
## memory than GSETS. Use of GSETS should be deprecated in the future.
message("[INIT] converting GSETS to list of integers...")
GSET_GENES <- sort(unique(unlist(GSETS)))  ## slow...
iGSETS <- parallel::mclapply(GSETS, function(a) match(a,GSET_GENES))  ## very slow!!!
names(iGSETS) <- names(GSETS)
## getGSETS <- function(gs) {lapply(iGSETS[gs],function(i) GSET_GENES[i]) }

## compare size savings!
object.size(GSETS) / 1e6
object.size(iGSETS) / 1e6

usethis::use_data(iGSETS)
usethis::use_data(GSET_GENES)
## usethis::use_data(getGSETS)   ## can a function be saved as data too???

##---------------------------------------------------------
message("[INIT] parsing collections...")
COLLECTIONS <- pgx.getGeneSetCollections(names(GSETS), min.size=10, max.size=99999)
COLLECTIONS <- COLLECTIONS[order(names(COLLECTIONS))]
usethis::use_data(COLLECTIONS)

##-----------------------------------------------------------------------------
## TISSUE/REFERENCE data sets
##-----------------------------------------------------------------------------
load(file.path(FILES,"sig/rna_tissue.rda"),verbose=1)  ## TISSUE and TISSUE.grp
usethis::use_data(TISSUE)
usethis::use_data(TISSUE_grp)

##-----------------------------------------------------------------------------
## Immune cell markers
##-----------------------------------------------------------------------------

## Really needed??? Seems already done by Nick
#IMMPROT <- read.csv(file.path(FILES,"sig/ImmProt-signature.csv"),row.names=1)
#IMMPROT_MARKERS <- rownames(read.csv(file.path(FILES,"sig/immprot-signature1000.csv"),row.names=1))
#DICE_MARKERS <- rownames(read.csv(file.path(FILES,"sig/DICE-signature1000.csv"),row.names=1))
#LM22 <- read.csv(file.path(FILES,"sig/LM22.txt"),sep="\t",row.names=1)
#LM22_MARKERS <- rownames(LM22)
#usethis::use_data(IMMPROT)
#usethis::use_data(TISSUE.grp)

##-----------------------------------------------------------------------------
## Colors
##-----------------------------------------------------------------------------

COLORS = rep(c(ggsci::pal_npg("nrc", alpha = 0.7)(10),
  ggsci::pal_aaas("default", alpha = 0.7)(10),
  ggsci::pal_d3("category10", alpha = 0.7)(10)),99)
BLUERED <- colorRampPalette(c("royalblue3","grey90","indianred3"))
PURPLEYELLOW <- colorRampPalette(c("purple","purple4","black","yellow4","yellow"))

usethis::use_data(COLORS)
usethis::use_data(BLUERED)
usethis::use_data(PURPLEYELLOW)
