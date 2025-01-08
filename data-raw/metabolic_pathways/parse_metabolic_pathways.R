library(dplyr)

#------------------------------------------------------------------------
# parse Reactome pathways
#------------------------------------------------------------------------
##setwd("~/Playground/playdata")

# Read file
data <- readr::read_tsv("./data-raw/metabolic_pathways/ChEBI2Reactome.txt")

# Inspect the data to determine the correct number of columns
head(data)
colnames(data) <- c("chebi", "stId", "reactome_url", "reactome_name", "evidence", "species")

head(data)
length(table(data$stId))

# keep Homo sapiens for now
table(data$species)
data <- data[data$species == "Homo sapiens", ]
head(data)
dim(data)

# append reactome path is to pathways name
data$reactome_name <- paste0(data$reactome_name, " [",data$stId,"]")

reactome_gmt <- tapply( data$chebi,  data$reactome_name, c)
reactome_gmt[[10]]
reactome_gmt <- lapply( reactome_gmt, function(m) paste0("CHEBI:",m))

#------------------------------------------------------------------------
# parse Wikipathways from AnnotationHub
#------------------------------------------------------------------------
library(AnnotationHub)

ah <- AnnotationHub()
query(ah, "wikipathways")
query(ah, "Homo sapiens")
qr <- query(ah, c("wikipathways", "metabolites"))
qr <- query(ah, c("wikipathways", "Homo sapiens", "metabolites"))
qr <- query(ah, c("wikipathways", "Homo sapiens"))
qr <- query(ah, c("wikipathways", "Homo sapiens", "ChEBI"))
qr <- query(ah, c("wikipathways", "Homo sapiens", "Pathways"))
#qr <- query(ah, c("pathbank", "Homo sapiens"))
qr
names(qr)

data <- qr[['AH91805']]
colnames(data)
head(data,10)
tail(data)

# append reactome path is to pathways name
# keep only numerical characters from ChEBI_ID
data$ChEBI_ID <- gsub("[^0-9]", "", data$ChEBI_ID)
##data$wp_name <- paste(data$pathway_name, data$wpid, sep = "_")
data$wp_name <- paste0(data$pathway_name," [", data$wpid, "]")

wiki_gmt <- tapply( data$ChEBI_ID, data$wpid, c)
wiki_gmt <- lapply( wiki_gmt, function(m) paste0("CHEBI:",sub("CHEBI:","",m)))
length(wiki_gmt)

## add gene/proteins of same pathway
wp.gset <- grep("_WP", rownames(playdata::GSETxGENE),value=1)
wp.gmt <- playbase::mat2gmt(t(playdata::GSETxGENE[wp.gset,]))
head(names(wp.gmt))
wp.gmt <- lapply( wp.gmt, function(p) paste0("SYMBOL:",p))
wp.gmt[[1]]

## merge
wp.id <- sub("_.*","",gsub(".*_WP","WP",names(wp.gmt)))
ii <- intersect(names(wiki_gmt), wp.id)
wp.gmt <- wp.gmt[match(names(wiki_gmt), wp.id)]

# merge genesets
##wiki_gmt <- playbase::read.gmt(gmt.file = "./data-raw/metabolic_pathways/metabolic_pathways_wp.gmt")
wiki_gmt <- mapply( union, wiki_gmt, wp.gmt )
length(wiki_gmt)
wp_id <- names(wiki_gmt)
wp_name <- data$pathway_name[match(wp_id, data$wpid)]
names(wiki_gmt) <- paste0(wp_name, " [", wp_id,"]")

#------------------------------------------------------------------------
# Multi-omics Pathways: Pathbank 2.0 (former SMDB)  IKnov2024
#------------------------------------------------------------------------

pw <- data.table::fread("~/Playground/public-db/pathbank.org/pathbank_pathways.csv")
mx <- data.table::fread("~/Playground/public-db/pathbank.org/pathbank_all_metabolites.csv")
px <- data.table::fread("~/Playground/public-db/pathbank.org/pathbank_all_proteins.csv")

head(pw)
head(mx)
dim(mx)
head(px)
colSums(!is.na(mx) & mx!="")


MX <- tapply( mx[1:100,"ChEBI ID"], mx[1:100,"PathBank ID"], c)
PATHBANK_MX <- tapply( mx[,"ChEBI ID"], mx[,"PathBank ID"], c)
PATHBANK_PX <- tapply( px[,"Uniprot ID"], px[,"PathBank ID"], c)
length(PATHBANK_MX)
length(PATHBANK_PX)

hist(sapply(PATHBANK_MX, length), breaks=50)
hist(sapply(PATHBANK_PX, length), breaks=50)

sum(sapply(PATHBANK_MX, length)>=3)
sum(sapply(PATHBANK_PX, length)>=3)
PATHBANK_MX <- PATHBANK_MX[ sapply(PATHBANK_MX, length)>=3 ]
PATHBANK_PX <- PATHBANK_PX[ sapply(PATHBANK_PX, length)>=3 ]

pp <- intersect(names(PATHBANK_MX), names(PATHBANK_PX))
length(pp)
PATHBANK_MX <- PATHBANK_MX[pp]
PATHBANK_PX <- PATHBANK_PX[pp]
PATHBANK_MX <- lapply( PATHBANK_MX, function(s) paste0("CHEBI:",s))
PATHBANK_PX <- lapply( PATHBANK_PX, function(s) paste0("UNIPROT:",s))
pathbank_gmt <- mapply(union, PATHBANK_PX, PATHBANK_MX )
pathbank_gmt[1:2]
length(pathbank_gmt)

head(pw)[,1:3]
ii <- match(names(pathbank_gmt), pw[["SMPDB ID"]])
smp.id <- names(pathbank_gmt)
pw.id <- pw[['PW ID']][match(smp.id, pw[['SMPDB ID']])]
pw.title <- pw[["Name"]][ii]
table(duplicated(pw.title))
pw.title <- paste0(pw.title," [",smp.id,"]")
head(pw.title)
names(pathbank_gmt) <- pw.title

pwsize <- sapply(pathbank_gmt, function(s) sum(grepl("CHEBI",s)))
pathbank_gmt <- pathbank_gmt[order(-pwsize)]

## take out duplicates
pwname <- gsub(" \\[.*","",names(pathbank_gmt))
table(duplicated(pwname))
pathbank_gmt <- pathbank_gmt[!duplicated(pwname)]
length(pathbank_gmt)

## translate UNIPROT to SYMBOL
mat <- playbase::gmt2mat( pathbank_gmt )
head(rownames(mat))
tail(rownames(mat))
head(colnames(mat))
ii <- grep("UNIPROT",rownames(mat))
pp <- sub("UNIPROT:","",rownames(mat)[ii])
symbol <- playbase::convert_probetype("human", pp, "SYMBOL", from_id="UNIPROT")
table(is.na(symbol))
rownames(mat)[ii] <- paste0("SYMBOL:",symbol)
mat <- mat[rownames(mat) != "SYMBOL:NA",]
mat <- mat[rownames(mat) != "CHEBI:NA",]
sum(duplicated(rownames(mat)))
mat <- mat[order(-Matrix::rowSums(mat)),]
mat <- mat[!duplicated(rownames(mat)),]
pathbank_gmt <- playbase::mat2gmt(mat)

length(pathbank_gmt)
pathbank_gmt[[2]]
ratio <- sapply( pathbank_gmt, function(g) mean(grepl("CHEBI",g)))
summary(ratio)


#------------------------------------------------------------------------
# combine pathways
#------------------------------------------------------------------------

pathways <- c(reactome_gmt, wiki_gmt, pathbank_gmt)
length(pathways)
head(names(pathways))

gmt1 <- list(CUSTOM = pathways)
pathways <- playbase::clean_gmt(gmt1, "METABOLOMICS")
head(names(pathways))

##pwsize <- sapply(pathways, length)
pwsize <- sapply(pathways, function(s) sum(grepl("CHEBI",s)))
pathways <- pathways[order(-pwsize)]
pwname <- gsub(" \\[.*","",names(pathways))
table(duplicated(pwname))
pathways <- pathways[!duplicated(pwname)]

pwsize <- sapply(pathways, length)
msize <- sapply(pathways, function(s) sum(grepl("CHEBI",s)))
table( msize >= 3)
table( msize >= 3 & pwsize>500)
pathways <- pathways[which(msize >=3)]

MSETxMETABOLITE <- playbase::createSparseGenesetMatrix(
    pathways,
    filter_genes = FALSE,
    min.geneset.size = 3,
    max.geneset.size = 500
)
dim(MSETxMETABOLITE)
head(rownames(MSETxMETABOLITE))
head(colnames(MSETxMETABOLITE))
tail(colnames(MSETxMETABOLITE))
colnames(MSETxMETABOLITE)

usethis::use_data(MSETxMETABOLITE, overwrite = TRUE)


