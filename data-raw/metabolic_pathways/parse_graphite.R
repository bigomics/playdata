#------------------------------------------------------------------------
# Multi-omics pathways edges/links from Graphite
#------------------------------------------------------------------------

##BiocManager::install("graphite")
library(graphite)
library(igraph)

  
if(0) {

  pathwayDatabases()
  pathwayDatabase()
  
  pathways <- graphite::pathways("hsapiens", "reactome")
  pathways <- graphite::pathways("hsapiens", "pathbank")
  pathways <- graphite::pathways("hsapiens", "smpdb")
  pathways <- graphite::pathways("hsapiens", "pharmgkb")
  length(pathways)
  names(pathways)[1:10]
  
  p <- pathways[[2]]
  p
  pathwayTitle(p)
  nodes(p)
  nodes(p, which = "mixed")
  edges(p)
  nodes(p, which = "mixed")
  edges(p, which = "mixed")
  head(edges(p, which = "metabolites"))
    
  g <- pathwayGraph(p, which = "mixed")
  print(g)
  g <- graph_from_graphnel(g)
  plot(g, vertex.label.cex=1, vertex.size=1)
  
}

pathwayDatabases()

## Collect all edges from all databases
all.db <- pathwayDatabases()[pathwayDatabases()$species=="mmusculus",2]
all.db <- pathwayDatabases()[pathwayDatabases()$species=="hsapiens",2]
all.db

all.db <- c("panther","pathbank","reactome","wikipathways")
db="reactome"
db="kegg"
all.edges <- c()

for(db in all.db) {
  pathways <- graphite::pathways("hsapiens", db)
  length(pathways)
  ##p1 = pathways[[k]]
  message("****** ",db," *******")
  
  getEE <- function(p1) {
    p1 <- convertIdentifiers(p1, "CHEBI")
    p1 <- convertIdentifiers(p1, "SYMBOL")
    ee <- edges(p1, which = "mixed")
    g <- pathwayGraph(p1, which = "mixed")
    g1 <- graph_from_graphnel(g)
    pw.id <- paste0(p1@database,":",p1@id)
    ee1 <- cbind( igraph::get.edgelist(g1), pw.id )
    ee1
  }

  message("length.pathways = ", length(pathways))
  ee <- parallel::mclapply( pathways, getEE, mc.cores = 48)
  ee <- do.call( rbind, ee )
  all.edges <- rbind(all.edges, ee)

  message("nrow.ee = ", nrow(all.edges))
  message("ncol.ee = ", ncol(all.edges))
  
  ## collapse multiple, count edges
  gr <- igraph::graph_from_edgelist(all.edges[,1:2], directed=FALSE)
  E(gr)$weight <- count_multiple(gr)
  sum(E(gr)$weight[!which_multiple(gr)])
  
  gr2 <- simplify(gr, edge.attr.comb = list(weight = "min"))
  ee <- igraph::get.edgelist(gr2) 
  ee <- apply(ee, 2, function(a) sub("SYMBOL:","",a))
  ee <- as.data.frame(ee) 
  ee$cost <- 1 / E(gr2)$weight
  GRAPHITE_PPI <- data.frame(from = ee[,1], to = ee[,2], cost = ee$cost)
  head(GRAPHITE_PPI)

  ## save(GRAPHITE_PPI, file="GRAPHITE_PPI.rda")
  ## load("~/Playground/public-db/pathbank.org/GRAPHITE_PPI.rda",verbose=TRUE)

  usethis::use_data(GRAPHITE_PPI, overwrite = TRUE)
  
  
}
