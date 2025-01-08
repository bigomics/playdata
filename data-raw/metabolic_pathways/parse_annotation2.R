library(dplyr)

##-----------------------------------------------------------------
# Create ID mapping table from different databases
##-----------------------------------------------------------------

# Read file (from https://ftp.ebi.ac.uk/pub/databases/chebi/). This
# database has ChEBI ID, name and definition of the metabolite. We
# restrict our annotation to only those molecules referred in ChEBI.
#setwd("~/Playground/playdata")
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
map <- metaboliteIDmapping::metabolitesMapping
object.size(map)/1e9
colnames(map)
colSums(!is.na(map))
map <- map[,c("ChEBI","HMDB","KEGG","CID","Name","CAS")]
colnames(map) <- sub("CID","PubChem",colnames(map))
tail(map)
dim(map)

# Notice not all ChEBI have mapping.
table(chebi$ID %in% map$ChEBI)
table(map$ChEBI %in% chebi$ID)

## match but not NA
match2 <- function(a,b)  ifelse(is.na(a), NA, match(a,b))

# add pathbank ID. We need this for the Pathbank SVG pathway images.
#system("cd ~/Downloads && wget https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip")
mx <- data.table::fread("~/Downloads/pathbank_all_metabolites.csv.zip")
sum(setdiff(map$HMDB,c(NA,"","-"))  %in% mx[["HMDB ID"]])
sum(setdiff(map$ChEBI,c(NA,"","-")) %in% mx[["ChEBI ID"]])
sum(setdiff(map$CAS,c(NA,"","-"))  %in% mx[["CAS"]])
i1 <- match2(map$HMDB, mx[["HMDB ID"]])
i1[is.na(map$HMDB)] <- NA
table(!is.na(i1))
i2 <- match2(map$ChEBI, mx[["ChEBI ID"]])
i2[is.na(map$ChEBI)] <- NA
table(!is.na(i2))
sum(!is.na(i2) & is.na(i1))
ii <- ifelse(!is.na(i2), i2, i1)
map$PATHBANK <- mx[["Metabolite ID"]][ii]

# Add REFMET. Get from  https://www.metabolomicsworkbench.org/databases/refmet/browse.php. It has nice annotation of super/main/sub class. Also has lipidmaps ID.
REFMET <- read.csv("./data-raw/metabolic_pathways/refmet_241218.csv")
dim(REFMET)
head(REFMET)
tail(sort(table(REFMET$refmet_name)))
table(REFMET$lipidmaps_id!="")

i1 <- match2(map$ChEBI, REFMET[["chebi_id"]])  ## best overlap
i2 <- match2(map$PubChem, REFMET[["pubchem_id"]])  ## best overlap
i3 <- match2(map$HMDB, REFMET[["hmdb_id"]])  ## best overlap
ii <- ifelse(!is.na(i1), i1, i2)
ii <- ifelse(!is.na(ii), ii, i3)
table(!is.na(ii))
map$refmet_name <- REFMET[["refmet_name"]][ii]
tail(sort(table(map$refmet_name)))

map$REFMET <- REFMET[["refmet_id"]][ii]
map$LIPIDMAPS <- REFMET[["lipidmaps_id"]][ii]

# convert all columns to ASCII character
map[] <- lapply(map, as.character)
map$Name <- trimws(iconv(map$Name,to="ascii//TRANSLIT"))
map$refmet_name <- trimws(iconv(map$refmet_name,to="ascii//TRANSLIT"))

# merge the two dataframes
ANNOTATION <- merge(chebi, map, by.x = "ID", by.y = "ChEBI", all.x = TRUE)
dim(ANNOTATION)
head(ANNOTATION)

# Add back HMDB entries with no ChEBI ID
ii <- which(!is.na(map$HMDB) & !map$HMDB %in% ANNOTATION$HMDB )
length(ii)
ANNOTATION <- merge(ANNOTATION, map[ii,], all.x = TRUE, all.y = TRUE)
dim(ANNOTATION)
rownames(ANNOTATION) <- NULL

# update missing CHEBI
ANNOTATION$ChEBI <- ifelse(!is.na(ANNOTATION$ChEBI), ANNOTATION$ChEBI, ANNOTATION$ID)
ANNOTATION$ID  <- ifelse(!is.na(ANNOTATION$ID), ANNOTATION$ID, ANNOTATION$ChEBI)

# remove unnecessary columns from annotation. This part only for ID
#conversion
cols <- c("ID","HMDB","PubChem","ChEBI","KEGG","PATHBANK","REFMET",
  "LIPIDMAPS","name","refmet_name")
cols %in% colnames(ANNOTATION)
ANNOTATION <- ANNOTATION[,cols]
head(ANNOTATION)
colSums(!is.na(ANNOTATION))
dim(ANNOTATION)

# Merge names. Notice we use AnnotHub name as default
ANNOTATION$name[ANNOTATION$name==""] <- NA
ANNOTATION$refmet_name[ANNOTATION$refmet_name==""] <- NA

ANNOTATION$name <- ifelse(!is.na(ANNOTATION$name),ANNOTATION$name, ANNOTATION$refmet_name)
ANNOTATION$refmet_name <- NULL

# remove those with no name???
no.name <- is.na(ANNOTATION$name) | ANNOTATION$name==""
table(no.name)
ANNOTATION <- ANNOTATION[!no.name,,drop=FALSE]  ## remove???

# remove duplicates
ANNOTATION <- ANNOTATION[order(rowSums(is.na(ANNOTATION))),]
dup1 <- (duplicated(ANNOTATION$ChEBI) & !is.na(ANNOTATION$ChEBI))
dup2 <- (duplicated(ANNOTATION$HMDB) & !is.na(ANNOTATION$HMDB))
sum(dup1)
sum(dup2)
table(!dup1 & !dup2)
ANNOTATION <- ANNOTATION[ (!dup1 & !dup2),]
dim(ANNOTATION)

# as rownames ChEBI or HMDB
table(is.na(ANNOTATION$ID))
annot.id <- ifelse(!is.na(ANNOTATION$ID), ANNOTATION$ID, ANNOTATION$HMDB)
table(is.na(annot.id))
ANNOTATION$ID <- annot.id
rownames(ANNOTATION) <- annot.id
ANNOTATION <- ANNOTATION[order(rownames(ANNOTATION)),]

head(ANNOTATION)
object.size(ANNOTATION)/1e9

## replace special-empty with NA
ANNOTATION <- apply(
  ANNOTATION,2, function(s) {sel=which(s %in% c('','-'));s[sel]=NA;s})
ANNOTATION <- data.frame(ANNOTATION)

## rename and save
METABOLITE_ID <- data.frame(ANNOTATION)
dim(METABOLITE_ID)
METABOLITE_ID$name <- NULL ##???
usethis::use_data(METABOLITE_ID, overwrite = TRUE)

##-----------------------------------------------------------------
# prepare metabolite metadata. This dataframe contains more
# information about metabolites like family, mass, formula, etc
##-----------------------------------------------------------------

metadata  <- ANNOTATION[,c("ID","name")]
metadata$definition <- chebi[match(ANNOTATION$ChEBI,chebi$ID),c("definition")]

i1 <- match2(metadata$ID, REFMET$chebi_id)
i2 <- match2(metadata$ID, REFMET$hmdb_id)
ii <- ifelse(is.na(i1),i2,i1)
table(is.na(ii))
metadata2 <- REFMET[ii,]

cols <- c('chebi_id','refmet_id','refmet_name','super_class','main_class','sub_class','formula','exactmass','pubchem_cid','hmdb_id','lipidmaps_id','kegg_id','inchi_key')
metadata2 <- metadata2[,cols]

METABOLITE_METADATA <- cbind(metadata, metadata2)
rownames(METABOLITE_METADATA) <- METABOLITE_METADATA$ID

dim(METABOLITE_METADATA)
head(METABOLITE_METADATA)
object.size(METABOLITE_METADATA) / 1e9
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
dim(METABOLITE_METADATA)
usethis::use_data(METABOLITE_METADATA, overwrite = TRUE)

head(METABOLITE_ANNOTATION)
head(METABOLITE_METADATA)
tail(METABOLITE_METADATA)

