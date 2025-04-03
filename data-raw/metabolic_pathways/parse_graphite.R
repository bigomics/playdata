#------------------------------------------------------------------------
# Multi-omics pathways edges/links from Graphite
#------------------------------------------------------------------------

##BiocManager::install(c("graphite","igraph"))
library(graphite)
library(igraph)

## Collect all edges from all databases: kegg, panther, pathbank,
## pharmgkb, reactome, smpdb, wikipathways
graphite::pathwayDatabases()
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


## convert identifiers to SYMBOL and CHEBI
load(file="graphite-pathways.rda", verbose=TRUE)
graphite::nodes(pathways[[1]][[1]], which='mixed')
names(pathways)
k="kegg"
for(k in names(pathways)) {
  pathways[[k]] <- lapply(pathways[[k]], convertIdentifiers, to="SYMBOL")
  pathways[[k]] <- lapply(pathways[[k]], convertIdentifiers, to="CHEBI")
}
save(pathways, file="graphite-pathways-symbol-chebi.rda")


## get all nodes to create gene sets. Get all edges for creating PPI.
load(file="graphite-pathways-symbol-chebi.rda",verbose=1)
names(pathways)
pathways[[1]][[1]]
all.nodes <- list()
all.edges <- list()
i="kegg"
for(i in names(pathways)) {
  nn <- lapply(pathways[[i]], function(p) graphite::nodes(p, which="mixed"))
  ee <- lapply(pathways[[i]], function(p) {
    e1 <- graphite::edges(p,which="mixed")
    cbind( paste0(e1$src_type,":",e1$src), paste0(e1$dest_type,":",e1$dest) )
  })  
  title <- sapply(pathways[[i]], function(p) p@title)
  id <- sapply(pathways[[i]], function(p) p@id)
  id <- sub("[:]","",id)
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
tail(names(pathways))
length(pathways)

##pwsize <- sapply(pathways, length)
wsize <- sapply(pathways, length)
msize <- sapply(pathways, function(s) sum(grepl("CHEBI",s)))
psize <- sapply(pathways, function(s) sum(grepl("SYMBOL",s)))
table( msize >= 3)
table( psize >= 200)
table( psize < 10)
table( msize >= 3 & psize >= 10)

sel <- which(msize >= 3 & psize >= 10 )
pathways <- pathways[sel]
length(pathways)

MSETxMETABOLITE <- playbase::createSparseGenesetMatrix(
    pathways,
    filter_genes = FALSE,
    min.geneset.size = 3,
    max.geneset.size = 1000,
    min_gene_frequency = 1
)

## There are many duplicated pathways with different names. We take
## out the duplicates and retain those with longer (more specific)
## name.
len.title <- nchar(rownames(MSETxMETABOLITE))
ii <- order(-Matrix::rowSums(MSETxMETABOLITE!=0),-len.title)
MSETxMETABOLITE <- MSETxMETABOLITE[ii,]
dim(MSETxMETABOLITE)

checksum <- apply(1*(MSETxMETABOLITE!=0),1,paste,collapse="")
table(duplicated(checksum))
names(checksum)[which(duplicated(checksum))]

short.title <- gsub("(TG|PC|CL|PE)\\(.*| \\[.*","",names(checksum))
table(duplicated(short.title))
  
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

##------------------------------------------------------------------------
## Create PPI edgelist. collapse multiple, count edges
##------------------------------------------------------------------------
library(igraph)

head(names(all.edges),1000)
tail(names(all.edges),1000)

## remove redundant (filtered out) pathways 
M <- MSETxMETABOLITE
sel.pathways <- unique(gsub("METABOLITE:","",rownames(M)))
sel.pathways <- sub("\\[hsa","[hsa:",sel.pathways)
sel <- which(names(all.edges) %in% sel.pathways )
length(sel)
sel.edges <- all.edges[sel]
length(sel.pathways)
length(all.edges)
length(sel.edges)

ee <- do.call(rbind, sel.edges)
sum(duplicated(ee))
gr <- igraph::graph_from_edgelist(ee, directed=FALSE)
E(gr)$weight <- count_multiple(gr)
sum(E(gr)$weight[!which_multiple(gr)])

hist(log10(1+E(gr)$weight), breaks=100)

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
## Extend metabolites adding first neighbours of proteins/metabolites
## in a geneset
##------------------------------------------------------------------------
source("~/Playground/playbase/dev/include.R",chdir=TRUE)

M <- MSETxMETABOLITE
G <- playdata::GSETxGENE
colnames(G) <- paste0("SYMBOL:",colnames(G))
head(colnames(G))
dim(G)

ppi <- GRAPHITE_PPI
extG  <- extend_metabolite_sets2(G, ppi, postfix="(ext2)", maxcost=0.1)
extM1 <- extend_metabolite_sets(M, ppi, postfix="(ext1)", add=FALSE, maxcost=0.1)
extM2 <- extend_metabolite_sets2(M, ppi, postfix="(ext2)", add=FALSE, maxcost=0.1)
M0 <- M
rownames(M0) <- paste(rownames(M0),"(ext0)")
extM <- Matrix::t(merge_sparse_matrix(Matrix::t(extG), Matrix::t(M0)))
extM <- Matrix::t(merge_sparse_matrix(Matrix::t(extM), Matrix::t(extM1)))
extM <- Matrix::t(merge_sparse_matrix(Matrix::t(extM), Matrix::t(extM2)))
dim(extM)
head(colnames(extM),100)
tail(colnames(extM))

## filter on size: minimum 3 metabolites, minimum 10 proteins
i1 <- grep("CHEBI:",colnames(extM))
i2 <- grep("SYMBOL:",colnames(extM))
num.mx <- Matrix::rowSums(extM[,i1]!=0)
num.px <- Matrix::rowSums(extM[,i2]!=0)
sel <- which(num.mx >= 3 & num.px >= 10)
extM <- extM[sel,]
dim(extM)

## write data object
XSETxMETABOLITE <- extM
usethis::use_data(XSETxMETABOLITE, overwrite = TRUE)

##load(file="../../data/XSETxMETABOLITE.rda",verbose=TRUE)

##------------------------------------------------------------------------
## Create igraphs
##------------------------------------------------------------------------
load(file="graphite-pathways-symbol-chebi.rda",verbose=1)
graphs <- list()
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
##------------------------------------------------------------------------
##------------------------------------------------------------------------

