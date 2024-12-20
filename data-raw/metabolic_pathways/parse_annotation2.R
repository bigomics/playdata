##setwd("~/Playground/playdata")
library(dplyr)

##-----------------------------------------------------------------
# Create ID mapping table from different databases
##-----------------------------------------------------------------

# Read file (from https://ftp.ebi.ac.uk/pub/databases/chebi/). This
# database has ChEBI ID, name and definition of the metabolite. We
# restrict our annotation to only those molecules referred in ChEBI.
chebi <- readr::read_tsv("./data-raw/metabolic_pathways/chebi_compounds_20240801_0501.tsv")
head(data.frame(chebi))
colnames(chebi)
chebi <- chebi[c("ID", "CHEBI_ACCESSION", "NAME", "DEFINITION")]
colnames(chebi) <- c("ID", "CHEBI_ACCESSION", "name", "definition")
head(chebi)
dim(chebi)
chebi <- apply(chebi, 2, function(x) trimws(iconv(x,to="ascii//TRANSLIT")))
chebi <- apply(chebi, 2, function(x) sub("null",NA,x))
chebi <- data.frame(chebi)

## Previously we used a files from metaboanalysts but need to
## download. Here we use the metaboliteIDmapping R package which also
## has more entries. 
library(metaboliteIDmapping)  # uses AnnotHub
map <- metabolitesMapping
colnames(map)
map <- map[,c("ChEBI","HMDB","KEGG","CID","Name")]
colnames(map)[4] <- "PubChem"
tail(map)
dim(map)

# Notice not all ChEBI have mapping.
table(chebi$ID %in% map$ChEBI)

# add pathbank ID. We need this for the Pathbank SVG pathway images.
#system("cd ~/Downloads && wget https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip")
mx <- data.table::fread("~/Downloads/pathbank_all_metabolites.csv.zip")
ii <- match(map$HMDB, mx[["HMDB ID"]])
table(!is.na(ii))
map$PATHBANK <- mx[["Metabolite ID"]][ii]

# Add REFMET. get from  https://www.metabolomicsworkbench.org/databases/refmet/browse.php
REFMET <- read.csv("./data-raw/metabolic_pathways/refmet_241218.csv")
dim(REFMET)
head(REFMET)
table(REFMET$lipidmaps_id!="")

##ii <- match(map$ChEBI, REFMET[["chebi_id"]])
length(intersect(map$PubChem, REFMET$pubchem_cid))
length(intersect(map$ChEBI, REFMET$chebi_id))
length(intersect(map$HMDB, REFMET$hmdb_id))
ii <- match(map$ChEBI, REFMET[["chebi_id"]])  ## best overlap
table(!is.na(ii))
map$refmet_name <- REFMET[["refmet_name"]][ii]
map$REFMET <- REFMET[["refmet_id"]][ii]
map$LIPIDMAPS <- REFMET[["lipidmaps_id"]][ii]

# convert all columns to character
map[] <- lapply(map, as.character)
sum(rowMeans(is.na(map))==1)
##map <- map[rowMeans(is.na(map[, -1])) < 1, ]
map$Name <- iconv(map$Name,to="ascii//TRANSLIT")

# merge the two dataframes
colnames(map)
colnames(chebi)
ANNOTATION <- merge(chebi, map, by.x = "ID", by.y = "ChEBI", all.x = TRUE)
dim(ANNOTATION)

# remove unnecessary columns from annotation. This part only for ID
#conversion
cols <- c("ID","HMDB","PubChem","CHEBI_ACCESSION","KEGG","PATHBANK","REFMET",
          "LIPIDMAPS","name","refmet_name")
cols %in% colnames(ANNOTATION)
ANNOTATION <- ANNOTATION[,cols]
head(ANNOTATION,20)
apply(ANNOTATION,2,function(s) length(unique(s)))

colnames(ANNOTATION) <- sub("CHEBI_ACCESSION","ChEBI",colnames(ANNOTATION))
ANNOTATION$ChEBI <- sub("CHEBI:","",ANNOTATION$ChEBI)

# remove those with no annotation
table(!is.na(ANNOTATION$name))
table(!is.na(ANNOTATION$refmet_name))
tail(sort(table(ANNOTATION$name)))
no.name <- is.na(ANNOTATION$name) | ANNOTATION$name==""
table(no.name)
ANNOTATION <- ANNOTATION[!no.name,,drop=FALSE]  ## remove???

# remove duplicates
ANNOTATION <- ANNOTATION[order(rowSums(is.na(ANNOTATION))),]
sum(duplicated(ANNOTATION$ChEBI))
ANNOTATION <- ANNOTATION[!duplicated(ANNOTATION$ChEBI),]
sum(duplicated(ANNOTATION$ID))
rownames(ANNOTATION) <- as.character(ANNOTATION$ChEBI)
ANNOTATION <- ANNOTATION[order(rownames(ANNOTATION)),]

##dim(playdata::METABOLITE_ANNOTATION) ## old
dim(ANNOTATION)
head(ANNOTATION)
object.size(ANNOTATION)
tail(sort(table(ANNOTATION$name)))

ANNOTATION <- apply(
  ANNOTATION,2, function(s) {sel=which(s %in% c('','-'));s[sel]=NA;s})
ANNOTATION <- data.frame(ANNOTATION)

ANNOTATION$name <- iconv(ANNOTATION$name,to="ascii//TRANSLIT")
ANNOTATION$refmet_name <- iconv(ANNOTATION$refmet_name,to="ascii//TRANSLIT")

METABOLITE_ANNOTATION <- data.frame(ANNOTATION)
dim(METABOLITE_ANNOTATION)
usethis::use_data(METABOLITE_ANNOTATION, overwrite = TRUE)

##-----------------------------------------------------------------
# prepare metabolite metadata. This dataframe contains more
# information about metabolites like family, mass, formula, etc
##-----------------------------------------------------------------

metadata1 <- chebi[match(ANNOTATION$ID,chebi$ID),]
metadata1$CHEBI_ACCESSION <- NULL

length(intersect(ANNOTATION$PubChem,REFMET$pubchem_cid))
length(intersect(ANNOTATION$ChEBI,REFMET$chebi_id))

chebi_id <- ANNOTATION$ChEBI
chebi_id[is.na(chebi_id)] <- "NA"
ii <- match(chebi_id, REFMET$chebi_id)
table(is.na(ii))
metadata2 <- REFMET[ii,]

cols <- c('chebi_id','refmet_id','refmet_name','super_class','main_class','sub_class','formula','exactmass','pubchem_cid','hmdb_id','lipidmaps_id','kegg_id','inchi_key')
metadata2 <- metadata2[,cols]
rownames(metadata2) <- rownames(ANNOTATION)

METABOLITE_METADATA <- cbind(metadata1, metadata2)
rownames(METABOLITE_METADATA) <- METABOLITE_METADATA$ID

dim(METABOLITE_METADATA)
head(METABOLITE_METADATA)
object.size(METABOLITE_METADATA)
tail(sort(table(METABOLITE_METADATA$name)))

## remove identifiers??
grep("_id",colnames(METABOLITE_METADATA),value=TRUE)
sel <- grep("_id|_cid",colnames(METABOLITE_METADATA))
METABOLITE_METADATA <- METABOLITE_METADATA[,-sel]
head(METABOLITE_METADATA)
colnames(METABOLITE_METADATA)

colSums(!is.na(METABOLITE_METADATA))
METABOLITE_METADATA <- apply(
  METABOLITE_METADATA,2, function(s) {sel=which(s %in% c('','-'));s[sel]=NA;s})
METABOLITE_METADATA <- data.frame(METABOLITE_METADATA)

METABOLITE_METADATA$name <- iconv(METABOLITE_METADATA$name,to="ascii//TRANSLIT")
METABOLITE_METADATA$refmet_name <- iconv(METABOLITE_METADATA$refmet_name,to="ascii//TRANSLIT")
METABOLITE_METADATA$definition <- iconv(METABOLITE_METADATA$definition,to="ascii//TRANSLIT")

colSums(!is.na(METABOLITE_METADATA))

usethis::use_data(METABOLITE_METADATA, overwrite = TRUE)

sort(table(METABOLITE_METADATA$super_class))
sort(table(METABOLITE_METADATA$main_class))
sort(table(METABOLITE_METADATA$sub_class))
