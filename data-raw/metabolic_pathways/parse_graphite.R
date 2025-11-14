##------------------------------------------------------------------------
## Multi-omics pathways edges/links from Graphite
##------------------------------------------------------------------------

## Creates:
##   MSETxMETABOLITE  - combined metaboliteset x metabolite sparsemat
##   GRAPHITE_PPI     - mixed PPI/MMI

##BiocManager::install(c("graphite","igraph"))
library(graphite)
library(igraph)
options(Ncpus=64)

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

## get all nodes to create gene sets. Get all edges for creating PPI.
load(file="graphite-pathways.rda",verbose=1)
names(pathways)
sapply(pathways, length)
pathways <- pathways[order(sapply(pathways, length))]
sapply(pathways, length)

all.nodes <- list()
all.edges <- list()
db="kegg"
db="pathpank"
db="wikipathways"
for(db in names(pathways)) {

  cat(">>> converting DB:", db, "<<<\n")
  pw <- pathways[[db]]
  #pw <- head(pw,10)
  #pw <- pw[c(1:5,436:438)]
  
  ##----------------------------------------------  
  ## Get all nodes as SYMBOL, CHEBI or LIPIDMAPS
  ##----------------------------------------------
  p0 <- graphite::convertIdentifiers(pw, to="SYMBOL")
  p1 <- graphite::convertIdentifiers(p0, to="CHEBI")
  p2 <- graphite::convertIdentifiers(p0, to="LIPIDMAPS")
  n0 <- lapply(p0, function(p) graphite::nodes(p))
  n1 <- lapply(p1, function(p) graphite::nodes(p, which='mixed'))
  n2 <- lapply(p2, function(p) graphite::nodes(p, which='mixed'))
  nodes <- mapply(union, n0, n1)
  nodes <- mapply(union, nodes, n2)

  ##----------------------------------------------
  ## get edges: proteins, metabolites, lipids
  ##----------------------------------------------  
  e0 <- lapply(p0, function(p) {
    e <- graphite::edges(p)
    cbind( paste0(e$src_type,":",e$src), paste0(e$dest_type,":",e$dest) )
  })  
  e1 <- lapply(p1, function(p) {
    e <- graphite::edges(p, which="mixed")
    cbind( paste0(e$src_type,":",e$src), paste0(e$dest_type,":",e$dest) )
  })  
  e2 <- lapply(p2, function(p) {
    e <- graphite::edges(p, which="mixed")
    cbind( paste0(e$src_type,":",e$src), paste0(e$dest_type,":",e$dest) )
  })  
  edges <- mapply(rbind, e0, e1)
  edges <- mapply(rbind, edges, e2)
  edges <- lapply( edges, function(e) e[!duplicated(e),,drop=FALSE] )

  ##----------------------------------------------
  ## extract title
  ##----------------------------------------------
  title <- sapply(pw, function(p) p@title)
  id <- sapply(pw, function(p) p@id)
  id <- sub("[:]","",id)
  db <- sapply(pw, function(p) p@database)    
  tt <- paste0(title, " [",id,"]")
  head(tt)
  names(nodes) <- tt
  names(edges) <- tt

  all.nodes <- c(all.nodes, nodes)
  all.edges <- c(all.edges, edges)

  save(all.nodes, file="graphite-nodes.rda")
  save(all.edges, file="graphite-edges.rda")
}

length(all.nodes)
head(names(all.nodes),20)

#------------------------------------------------------------------------
# Create pathways gene sets
#------------------------------------------------------------------------
load(file="graphite-nodes.rda", verbose=1)
load(file="graphite-edges.rda", verbose=1)
length(all.nodes)
tail(names(all.nodes))
tail(all.nodes)

ngene  <- sapply(all.nodes, function(n) sum(grepl("SYMBOL",n)))
nlipid <- sapply(all.nodes, function(n) sum(grepl("LIPID",n)))
nchebi <- sapply(all.nodes, function(n) sum(grepl("CHEBI",n)))

table(nchebi>=3)
table(nlipid>=3)
table(ngene>=10)
table(ngene>=3)
table(ngene>=3 & nchebi>=3)
table(has.lipid, has.chebi)

## cleanup names
gmt <- list(CUSTOM = all.nodes)
gmt <- playbase::clean_gmt(gmt, "METABOLITE")
head(names(gmt))
tail(names(gmt))
length(gmt)

##pwsize <- sapply(gmt, length)
wsize <- sapply(gmt, length)
msize <- sapply(gmt, function(s) sum(grepl("CHEBI",s)))
lsize <- sapply(gmt, function(s) sum(grepl("LIPID",s)))
psize <- sapply(gmt, function(s) sum(grepl("SYMBOL",s)))
sel <- ((msize >= 3 | lsize >= 3) & psize >= 3)
table(sel)
gmt <- gmt[sel]
length(gmt)

MSETxMETABOLITE <- playbase::createSparseGenesetMatrix(
    gmt,
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

## we use a checksum to see if rows are duplicated
checksum <- apply(1*(MSETxMETABOLITE!=0),1,paste,collapse="")
table(duplicated(checksum))

## Check if short titles are duplicated (??). 
short.title <- gsub("(TG|PC|CL|PE)\\(.*","",names(checksum))
table(duplicated(short.title))

## The SMPDB database has many duplicated pathways with different
## names. This is very annoying.
is.smp <- grepl("SMP",names(checksum))
table(is.smp & duplicated(checksum) & duplicated(short.title))
sel <- (!(is.smp & duplicated(checksum) & duplicated(short.title)))
table(sel)

MSETxMETABOLITE <- MSETxMETABOLITE[sel,]
dim(MSETxMETABOLITE)
dim(playdata::MSETxMETABOLITE)

## write data object
usethis::use_data(MSETxMETABOLITE, overwrite = TRUE)

##------------------------------------------------------------------------
## Create PPI edgelist. collapse multiple, count edges
##------------------------------------------------------------------------
library(igraph)

head(names(all.edges),100)
tail(names(all.edges),100)

## remove redundant (filtered out) pathways 
M <- MSETxMETABOLITE
dim(M)
sel.pathways <- unique(gsub("METABOLITE:","",rownames(M)))
sel.pathways <- gsub("[:]","",sel.pathways)
names.all.edges <- gsub("[:]","",names(all.edges))
sel <- which(names.all.edges %in% sel.pathways )
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
dim(GRAPHITE_PPI)
head(GRAPHITE_PPI)
usethis::use_data(GRAPHITE_PPI, overwrite = TRUE)

##------------------------------------------------------------------------
## Create igraphs
##------------------------------------------------------------------------
load(file="graphite-pathways.rda",verbose=1)
graphs <- list()
i=j=1
for(i in names(pathways)) {
  pw <- pathways[[i]]
  p0 <- graphite::convertIdentifiers(pw, to="SYMBOL")
  p1 <- graphite::convertIdentifiers(p0, to="CHEBI")
  p2 <- graphite::convertIdentifiers(p0, to="LIPIDMAPS")
  for(j in 1:length(pw)) {
    e <- pw@entries[[j]]
    tt <- paste0(e@database,":",e@title, " [", e@id, "]")
    g1 <- pathwayGraph(p1[[j]], which="mixed")
    g2 <- pathwayGraph(p2[[j]], which="mixed")
    ig1 <- igraph::graph_from_graphnel(g1)
    ig2 <- igraph::graph_from_graphnel(g2)
    graphs[[tt]] <- igraph::union(ig1,ig2)
  }
}
save(graphs, file="graphite-igraphs.rda")


##------------------------------------------------------------------------
##------------------------------------------------------------------------
##------------------------------------------------------------------------

