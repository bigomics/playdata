
##setwd("..")

## From https://epsd.biocuckoo.cn/Download/Total.zip
system("curl https://epsd.biocuckoo.cn/Download/Total.zip --output ~/Downloads/Total.zip")
P1 <- data.table::fread(unzip("~/Downloads/Total.zip","Total.txt"))
colnames(P1)
P1 <- P1[,c("UniProt ID","Position","AA")]
colnames(P1) <- c("UniProt","Position","Residue")
P1 <- as.data.frame(P1)
sum(duplicated(P1[,1:2]))
P1 <- P1[!duplicated(P1[,1:2]),]

## https://dbpsp.biocuckoo.cn/Download/PhosphorylationData/Total.txt
P2 <- data.table::fread("https://dbpsp.biocuckoo.cn/Download/PhosphorylationData/Total.txt")
P2 <- P2[,c("Uniprot.ID","Position","Residue")]
##P2 = read.csv("data-raw/extdata/phospho/phosphosites-dbpsp_20190701.csv.gz")
P2 <- as.data.frame(P2)

dim(P1)
dim(P2)
head(P1)
head(P2)

p2.names <- paste0(P2[,1],"_",P2[,2])

rownames(P1) <- paste0(P1[,1],"_",P1[,2])
rownames(P2) <- paste0(P2[,1],"_",P2[,2])
length(intersect(rownames(P2),rownames(P1)))

colnames(P2) <- colnames(P1)
PHOSPHOSITE <- rbind(P1, P2)
rownames(PHOSPHOSITE) <- NULL  ## saves memory!
dim(PHOSPHOSITE)
usethis::use_data(PHOSPHOSITE, overwrite = TRUE)
