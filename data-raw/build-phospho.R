
##setwd("..")

P1 = read.csv("data-raw/extdata/phospho/phosphosites-dbpaf_20170710.csv")
P2 = read.csv("data-raw/extdata/phospho/phosphosites-dbpsp_20190701.csv")

dim(P1)
dim(P2)
head(P1)
head(P2)
rownames(P1) <- paste0(P1[,1],"_",P1[,2])
rownames(P2) <- paste0(P2[,1],"_",P2[,2])
pp <- intersect(rownames(P2),rownames(P1))
length(pp)
head(P1[pp,],10)
head(P2[pp,],10)

P2 <- cbind(P2, Species = NA)
colnames(P2) <- colnames(P1)
PHOSPHOSITE <- rbind(P1, P2)
rownames(PHOSPHOSITE) <- NULL  ## saves memory!
usethis::use_data(PHOSPHOSITE, overwrite = TRUE)
