##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## NOTE: This was the code piece creating the global-init.rda or
## lib/sysdata.rda file that was loaded at app start. See still the
## loading of sysdata.rda in
## omicsplayground/components/app/R/global.R. This loaded all data in
## memory. If playdata works, loading sysdata.rda should be not necessary
## anymore.

## NOTE: Do not run this file for playdata. This is just a reference
## to see how the data objects were created. Ideally, we should make a
## new script that recreates the data/*rda objects from raw csv as
## described below.
##


##-----------------------------------------------------------------------------
## GLOBAL variables
##-----------------------------------------------------------------------------
if(0) {
    RDIR='~/Playground/omicsplayground/R'
    FILES='~/Playground/omicsplayground/lib'
    source(file.path(RDIR,"pgx-include.R"),local=TRUE)  ## pass local vars

    FILES = "~/Playground/omicsplayground/lib"
}
source(file.path(RDIR,"pgx-functions.R"),local=TRUE)  ## pass local vars

## Caching the init files
INIT.FILE <- file.path(FILES,"global-init.rda")
INIT.FILE <- "../cache/global-init.rda" ## avoid rw permission
##unlink(INIT.FILE)
INIT.FILE

file.exists(INIT.FILE)

if(1 && file.exists(INIT.FILE)) {    

    message("[INIT] loading cached INIT file ",INIT.FILE)
    t0 <- Sys.time()
    load(INIT.FILE, verbose=1)
    message("Loading cache took: ", round(Sys.time() - t0), " seconds")

} else {
    
    message("[INIT] no INIT file! building INIT from scratch.")
    message("[INIT] INIT.FILE = ", INIT.FILE)    
    t0 <- Sys.time()

    oldvars <- ls()

    ## All gene families in Human UPPER CASE
    require(org.Hs.eg.db)
    GENE.TITLE  = unlist(as.list(org.Hs.egGENENAME))
    GENE.SYMBOL = unlist(as.list(org.Hs.egSYMBOL))
    names(GENE.TITLE) = GENE.SYMBOL
    ##GSET.PREFIX.REGEX = paste(paste0("^",GSET.PREFIXES,"_"),collapse="|")
    GSET.PREFIX.REGEX="^BIOCARTA_|^C2_|^C3_|^C7_|^CHEA_|^GOBP_|^GOCC_|^GOMF_|^HALLMARK_|^KEA_|^KEGG_|^PID_|^REACTOME_|^ST_"
    GENE.SUMMARY = read.csv(file.path(FILES,"gene-summary.csv"),row.names=1)
    GENE.SUMMARY = array(GENE.SUMMARY[,1], dimnames=list(rownames(GENE.SUMMARY)))
    
    ## GENExGENE <- readRDS(file=file.path(FILES,"GENExGENE-cosSparseKNN500-XL.rds"))
    GSETxGENE <- readRDS(file.path(FILES,"gset-sparseG-XL.rds"))
    load(file.path(FILES,"gmt-all.rda"),verbose=1)
    GSETS = gmt.all;remove(gmt.all)

    message("[INIT] parsing gene families...")
    FAMILIES <- pgx.getGeneFamilies(GENE.SYMBOL, FILES=FILES, min.size=10, max.size=9999)
    fam.file <- file.path(FILES,"custom-families.gmt")
    if(file.exists(fam.file)) {
        custom.gmt = read.gmt(file.path(FILES,"custom-families.gmt"),add.source=TRUE)
        names(custom.gmt)
        FAMILIES= c(FAMILIES, custom.gmt)
    }
    FAMILIES[["<all>"]] <- GENE.SYMBOL
    f1 <- FAMILIES
    names(f1) <- paste0("FAMILY:",names(f1))
    names(f1) <- sub("FAMILY:<all>","<all>",names(f1))
    GSETS <- c(GSETS,f1)

    ## convert to integer list (more efficient)
    message("[INIT] converting GSETS to list of integers...")
    GSET.GENES <- sort(unique(unlist(GSETS)))  ## slow...
    iGSETS <- parallel::mclapply(GSETS, function(a) match(a,GSET.GENES))  ## slow...
    names(iGSETS) <- names(GSETS)
    getGSETS <- function(gs) {
        lapply(iGSETS[gs],function(i) GSET.GENES[i])
    }
        
    message("[INIT] parsing collections...")
    COLLECTIONS <- pgx.getGeneSetCollections(names(GSETS), min.size=10, max.size=99999)
    COLLECTIONS <- COLLECTIONS[order(names(COLLECTIONS))]

    remove(list=c("custom.gmt","f1","GSETS"))

    ##-----------------------------------------------------------------------------
    ## TISSUE/REFERENCE data sets
    ##-----------------------------------------------------------------------------
    load(file.path(FILES,"sig/rna_tissue.rda"))  ## TISSUE and TISSUE.grp

    ##-----------------------------------------------------------------------------
    ## Immune cell markers
    ##-----------------------------------------------------------------------------

    ## Really needed???
    IMMPROT <- read.csv(file.path(FILES,"sig/ImmProt-signature.csv"),row.names=1)
    IMMPROT_MARKERS <- rownames(read.csv(file.path(FILES,"sig/immprot-signature1000.csv"),row.names=1))
    DICE_MARKERS <- rownames(read.csv(file.path(FILES,"sig/DICE-signature1000.csv"),row.names=1))
    LM22 <- read.csv(file.path(FILES,"sig/LM22.txt"),sep="\t",row.names=1)
    LM22_MARKERS <- rownames(LM22)

    ##-----------------------------------------------------------------------------
    ## Meta MA-profiles (fold changes) of all experiments
    ##-----------------------------------------------------------------------------
    ##load( file.path(FILES,"allMA-pub.rda"), verbose=1)
    ##load(file.path(FILES,"allFoldChanges-pub-8k.rda"))
    ##PROFILES <- list(M=allM, A=allA, FC=allFC)
    ##allFC <- pgx.readDatasetProfiles(PGX.DIR, file="datasets-allFC.csv") 
    ##PROFILES <- list(FC=allFC)
    ##remove(allA)
    ##remove(allM)
    ##remove(allFC)

    ##-----------------------------------------------------------------------------
    ## Colors
    ##-----------------------------------------------------------------------------
        
    COLORS = rep(RColorBrewer::brewer.pal(8,"Set2"),99)
    COLORS = rep(c(ggsci::pal_npg("nrc", alpha = 0.7)(10),
                   ggsci::pal_aaas("default", alpha = 0.7)(10),
                   ggsci::pal_d3("category10", alpha = 0.7)(10)),99)
##    BLUERED <- grDevices::colorRampPalette(
##        rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#EEEEEE",
##              "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))
##    BLUERED <- function(n=64) suppressWarnings(gplots::colorpanel(n,low="royalblue3",mid="grey90",high="indianred3"))
    BLUERED <- colorRampPalette(c("royalblue3","grey90","indianred3"))
    PURPLEYELLOW <- colorRampPalette(c("purple","purple3","black","yellow3","yellow"))
    PURPLEYELLOW <- colorRampPalette(c("purple","purple4","black","yellow4","yellow"))

    newvars <- setdiff(ls(), oldvars)
    newvars

    message("Creating global init took: ", round(Sys.time() - t0), " seconds")
    message("[INIT] saving INIT file ", INIT.FILE)    
    save( list=newvars, file=INIT.FILE)
    
}
