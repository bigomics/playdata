##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

path_to_gmt <- "data-raw/extdata/gmt/"

##----------------------------------------------------------------------
##--------------- cell type signatures from ImSig ----------------------
##----------------------------------------------------------------------
S <- read.csv("data-raw/extdata/imsig-signatures-genes.csv", skip=2)

imsig.gmt <- playbase::convert.gmt(S$Gene.Symbol, S$ImSig.cell.type)
playbase::write.gmt(imsig.gmt, file= paste0(path_to_gmt,"celltype_imsig.gmt"))

##----------------------------------------------------------------------
##----------- cell type signatures from xCELL (collapsed) --------------
##----------------------------------------------------------------------

xcell.data <- playbase::read.gmt("data-raw/extdata/xCell_celltype_signatures.txt")
cell.type <- sub("_.*","",names(xcell.data))
table(cell.type)
xcell.gmt <- tapply( xcell.data, cell.type, function(s) sort(unique(unlist(s))))
playbase::write.gmt(xcell.gmt, file=paste0(path_to_gmt,"celltype_xcell.gmt"))

##----------------------------------------------------------------------
##----------- build the GMT-all object (all gene sets) -----------------
##----------------------------------------------------------------------

##tt = read.gmt("../../Data/Enrichr/Tissue_Jensen.txt")
##head(sapply(tt,length),100)

require(parallel)
gmt.files2 = dir(path_to_gmt, pattern=".gmt$|.txt$", full.names=TRUE)

# remove unwanted files

gsets_irrelevant <- c(
    "GTEx_Tissue_Expression",
    "Epigenomics_Roadmap_HM",
    "NIH_Funded_PIs",
    "Genes_Associated_with_NIH_Grants",
    "DrugMatrix",
    "MSigDB_Computational",
    "MSigDB_Hallmark",
    "MSigDB_Oncogenic",
    "HDSigDB",
    "ENCODE",
    "Old_CMAP",
    "Enrichr",
    "Rare",
    "LINCS",
    "C1",
    "C2",
    "C3",
    "C4",
    "C5",
    "C6",
    "C7",
    "C8",
    "Orphanet",
    "DrugMatrix",
    "SILAC_Phosphoproteomics",
    "Proteomics_Drug_Atlas",
    "ProteomicsDB",
    "CCLE_Proteomics",
    "SysMyo",
    "GeneSigDB",
    "CellMarker_Augmented",
    "PFOCR",
    "Elsevier_Pathway_Collection",
    "BioCarta_2016"
    )

pattern_irrelevant <- paste(gsets_irrelevant, collapse = "|")
gmt.files2 <- gmt.files2[grep(pattern_irrelevant, gmt.files2, invert=TRUE)]

gmt.all = mclapply(gmt.files2, playbase::read.gmt)

names(gmt.all) = gmt.files2
names(gmt.all) = gsub(".*/|.txt$|.gmt$", "", names(gmt.all))

# add drug prefix in certain datasets
gsets_drug_related <- c(
    "DSigDB",
    "RNA-Seq_Disease_Gene",
    "IDG_Drug_Targets")

pattern_to_add_prefix <- paste(gsets_drug_related, collapse = "|")

idx_to_modify <- grep(pattern_to_add_prefix, names(gmt.all))

names(gmt.all)[idx_to_modify] <- paste("Drug", names(gmt.all)[idx_to_modify], sep = "_")

# end of adding drug prefix 

# convert bioplex hugo symbol to gene id, in order to access database

ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Function to get the gene ID from HGNC
getGeneID <- function(gene){
    gene_info <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "entrezgene_id"),
        filters = "hgnc_symbol",
        values = gene,
        mart = ensembl)
        return(gene_info)
}

bioplex_ds_id <- grep("BioPlex",names(gmt.all))

bioplex_gene_ids <- data.frame(hgnc = names(gmt.all[[bioplex_ds_id]]))

bioplex_gene_mapping <- getGeneID(bioplex_gene_ids$hgnc)


bioplex_gene_ids$entrezgene_id <- bioplex_gene_mapping$entrezgene_id[match(bioplex_gene_ids$hgnc, bioplex_gene_mapping$hgnc_symbol)] 

bioplex_gene_ids$entrezgene_id <- as.character(bioplex_gene_ids$entrezgene_id)
bioplex_gene_ids$comb <- apply(bioplex_gene_ids, 1, paste,collapse = "_")
bioplex_gene_ids$comb <- sub("_NA","",bioplex_gene_ids$comb)

names(gmt.all[[bioplex_ds_id]]) <- bioplex_gene_ids$comb

# end of convert bioplex hugo symbol to gene id, in order to access database

gmt.db <- gsub("[_.-].*|.txt$|.gmt$", "", names(gmt.all))

gmt.all <- playbase::clean_gmt(gmt.all, gmt.db)

table(sub(":.*","",names(gmt.all)))

# remove C2:REACTOME as duplicate

gmt.all <-  gmt.all[grep("C2:REACTOME", names(gmt.all), invert = TRUE)]
length(gmt.all)
save(gmt.all, file="data-raw/extdata/gmt-all.rda")
saveRDS(gmt.all, file="data-raw/extdata/gmt-all.rds")


## NEED RETHINK!!!! this should be improved using BioMart...
cat("Converting gsets to mouse ID...\n")
mouse.genes = as.character(unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL)))
names(mouse.genes) = toupper(mouse.genes)
gmt.all <- mclapply(gmt.all[], function(s) setdiff(as.character(mouse.genes[s]),NA), mc.cores=1)
save(gmt.all, file="data-raw/extdata/gmt-all-mouse.rda")
saveRDS(gmt.all, file="data-raw/extdata/gmt-all-mouse.rds")
remove(gmt.all)