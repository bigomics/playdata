#------------------------------------------------------------------------
# Multi-omics pathways edges/links from Graphite
#------------------------------------------------------------------------

##BiocManager::install("graphite")
library(graphite)
library(igraph)

## Collect all edges from all databases: kegg, panther, pathbank,
## pharmgkb, reactome, smpdb, wikipathways
pathwayDatabases()
##dbs <- graphite:::.dbs[['hsapiens']]
dbs <- pathwayDatabases()[pathwayDatabases()$species=="hsapiens",2]
dbs
all.db <- setdiff(dbs, "smpdb")
pathways <- list()
db="kegg"
for(db in all.db) {
  paths <- graphite::pathways("hsapiens", db)
    pathways[[db]] <- paths
}
save(pathways, file="graphite-pathways.rda")

## get all nodes to create gene sets. Get all edges for creating PPI.
load(file="graphite-pathways-symbol-chebi.rda",verbose=1)
names(pathways)
pathways[[1]][[1]]
all.nodes <- list()
all.edges <- list()
for(i in names(pathways)) {
  nn <- lapply(pathways[[i]], function(p) graphite::nodes(p, which="mixed"))
  ee <- lapply(pathways[[i]], function(p) {
    e1 <- graphite::edges(p,which="mixed")
    cbind( paste0(e1$src_type,":",e1$src), paste0(e1$dest_type,":",e1$dest) )
  })  
  title <- sapply(pathways[[i]], function(p) p@title)
  id <- sapply(pathways[[i]], function(p) p@id)
  db <- sapply(pathways[[i]], function(p) p@database)    
  tt <- paste0(title, " [",id,"]")
  head(tt)
  names(nn) <- tt
  names(ee) <- tt
  all.nodes <- c(all.nodes, nn)
  all.edges <- c(all.edges, ee)
}
save(all.nodes, file="graphite-nodes.rda")
save(all.edges, file="graphite-edges.rda")

load(file="graphite-nodes.rda", verbose=1)
load(file="graphite-edges.rda", verbose=1)

#------------------------------------------------------------------------
# Create pathways gene sets
#------------------------------------------------------------------------

## cleanup names
gmt1 <- list(CUSTOM = all.nodes)
pathways <- playbase::clean_gmt(gmt1, "METABOLITE")
head(names(pathways))

##pwsize <- sapply(pathways, length)
msize <- sapply(pathways, function(s) sum(grepl("CHEBI",s)))
pathways <- pathways[order(-msize)]
pwname <- gsub(" \\[.*","",names(pathways))
table(duplicated(pwname))
pathways <- pathways[!duplicated(pwname)]

wsize <- sapply(pathways, length)
msize <- sapply(pathways, function(s) sum(grepl("CHEBI",s)))
psize <- sapply(pathways, function(s) sum(grepl("SYMBOL",s)))
table( msize >= 3)
pathways <- pathways[which(msize >=3)]
length(pathways)

MSETxMETABOLITE <- playbase::createSparseGenesetMatrix(
    pathways,
    filter_genes = FALSE,
    min.geneset.size = 3,
    max.geneset.size = 1000
)

## There are many duplicated pathways with different names. We take
## out the duplicates and retain those with longer (more specific)
## name.
len.title <- nchar(rownames(MSETxMETABOLITE))
ii <- order(-Matrix::rowSums(MSETxMETABOLITE!=0),-len.title)
MSETxMETABOLITE <- MSETxMETABOLITE[ii,]

dim(MSETxMETABOLITE)
checksum <- sapply(apply(MSETxMETABOLITE!=0,1,which),sum)
wdup <- which(duplicated(checksum))  
sum(duplicated(checksum))
sum(!duplicated(checksum))

short.title <- gsub("(TG|PC|CL|PE)\\(.*| \\[.*","",names(checksum))
sum(duplicated(short.title))
  
sel <- which(!duplicated(checksum) & !duplicated(short.title))
length(sel)
names(checksum)[sel]
MSETxMETABOLITE <- MSETxMETABOLITE[sel,]

dim(MSETxMETABOLITE)
dim(playdata::MSETxMETABOLITE)
head(rownames(MSETxMETABOLITE))
head(colnames(MSETxMETABOLITE))
tail(colnames(MSETxMETABOLITE))

## write data object
usethis::use_data(MSETxMETABOLITE, overwrite = TRUE)

head(colnames(playdata::MSETxMETABOLITE))
tail(colnames(playdata::MSETxMETABOLITE))

##------------------------------------------------------------------------
## Create PPI edgelist. collapse multiple, count edges
##------------------------------------------------------------------------
library(igraph)

ee <- do.call(rbind, all.edges)
sum(duplicated(ee))
gr <- igraph::graph_from_edgelist(ee, directed=FALSE)
E(gr)$weight <- count_multiple(gr)
sum(E(gr)$weight[!which_multiple(gr)])

gr2 <- simplify(gr, edge.attr.comb = list(weight = "min"))
ee <- igraph::as_edgelist(gr2) 
ee <- apply(ee, 2, function(a) sub("SYMBOL:","",a))
ee <- as.data.frame(ee) 
ee$cost <- 1 / E(gr2)$weight
GRAPHITE_PPI <- data.frame(from = ee[,1], to = ee[,2], cost = ee$cost)
head(GRAPHITE_PPI)

## save(GRAPHITE_PPI, file="GRAPHITE_PPI.rda")
## load("~/Playground/public-db/pathbank.org/GRAPHITE_PPI.rda",verbose=TRUE)
usethis::use_data(GRAPHITE_PPI, overwrite = TRUE)

##------------------------------------------------------------------------
## Create PPI edgelist. collapse multiple, count edges
##------------------------------------------------------------------------

igraphs <- list()
i=j=1
for(i in names(pathways)) {
  pw <- pathways[[i]]
  for(j in 1:length(pw)) {
    gr <- pathwayGraph(pw[[j]], which="mixed")
    e <- pw@entries[[j]]
    tt <- paste0(e@database,":",e@title, " [", e@id, "]")
    graphs[[tt]] <- igraph::graph_from_graphnel(gr)
  }
}
save(graphs, file="graphite-igraphs.rda")

##------------------------------------------------------------------------
## Extend genesets
##------------------------------------------------------------------------
## library(igraph)

## ## write data object
## M <- playdata::MSETxMETABOLITE
## dim(M)
## ppi <- playdata::GRAPHITE_PPI
## ppi[,1] <- ifelse( grepl("CHEBI",ppi[,1]), ppi[,1], paste0("SYMBOL:",ppi[,1]))
## ppi[,2] <- ifelse( grepl("CHEBI",ppi[,2]), ppi[,2], paste0("SYMBOL:",ppi[,2]))
## sel <- which( grepl("CHEBI",ppi[,1]) & grepl("CHEBI",ppi[,2]) & ppi[,3] <= 0.33 )
## gr <- graph_from_edgelist( as.matrix(ppi[sel,1:2]) )
## ppi_mat <- as.matrix(gr)
## dim(ppi_mat)
## table(colnames(M) %in% rownames(ppi_mat))
## ## build neigborhood matrix for metabolites
## ppi0 <- rbind("na"=0, cbind("na"=0, ppi_mat))
## ii <- match(colnames(M),rownames(ppi0))
## ii[is.na(ii)] <- 1
## B <- ppi0[ii,ii]
## colnames(B)=rownames(B)=colnames(M)
## diag(B) <- 1
## ## propagate neighbor
## M1 <- M %*% B
## write data object
#MSETxMETABOLITE2 <- M1
#usethis::use_data(MSETxMETABOLITE2, overwrite = TRUE)

