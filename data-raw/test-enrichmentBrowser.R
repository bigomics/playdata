# BiocManager::install("EnrichmentBrowser")
# BiocManager::install("msigdbr")
library(EnrichmentBrowser)


# showAvailableSpecies(db = c("go", "kegg", "msigdb", "enrichr"), cache = TRUE)
# showAvailableCollections(org, db = c("go", "kegg", "msigdb", "enrichr"), cache = TRUE)

db.list <- c("go", "msigdb", "enrichr")
db.species <- lapply(db.list, EnrichmentBrowser::showAvailableSpecies, cache=TRUE)

db.collections <- lapply(
  db.list,function(db){
    db_collection <- showAvailableCollections("hsa", db=db))
    write.csv(db.collection, file=paste0("data-raw//",db,".csv"))
    return(db_collection)
}

idTypes("hsa")

gs <- getGenesets(org = "hsa", gene.id.type = "SYMBOL", db = "enrichr", lib="WikiPathway_2021_Human")
gs[[1]]


gs <- getGenesets(org = "cel", gene.id.type = "SYMBOL", db = "go", onto="BP", mode="GO.db")
names(gs)


writeGMT(gs, gmt.file)



##https://rdrr.io/bioc/EnrichmentBrowser/man/getGenesets.html
go.gs <- getGenesets(org = "hsa", db = "go", onto = "BP", mode = "GO.db")
#go.gs <- getGenesets(org = "hsa", db = "go", mode = "biomart")

# (3) Obtaining *H*allmark gene sets from MSigDB
showAvailableSpecies(db = "msigdb")
showAvailableCollections(db = "msigdb") 
hall.gs <- getGenesets(org = "hsa", db = "msigdb", cat = "H")

# list supported species for obtaining gene sets from MSigDB

# (4) Obtaining gene sets from Enrichr
tfppi.gs <- getGenesets(org = "hsa", db = "enrichr", lib = "Transcription_Factor_PPIs")

# list supported species for obtaining gene sets from Enrichr
showAvailableSpecies(db = "enrichr")
showAvailableCollections(org = "hsa", db = "enrichr")        

# (6) parsing gene sets from GMT
gmt.file <- system.file("extdata/hsa_kegg_gs.gmt",package = "EnrichmentBrowser")
gs <- getGenesets(gmt.file)     

# (7) writing gene sets to file
writeGMT(gs, gmt.file)



##---------------------------------------------------------------------------
##---------------------------------------------------------------------------
##---------------------------------------------------------------------------


## (1) simulating expression data: 100 genes, 12 samples
se <- makeExampleData(what="SE") 
exprsHeatmap(expr=assay(se), grp=as.factor(se$GROUP))
se <- deAna(se)
pdistr(rowData(se)$ADJ.PVAL)
volcano(fc=rowData(se)$FC, p=rowData(se)$ADJ.PVAL)


  # (2) gene sets:
# draw 10 gene sets with 15-25 genes
gs <- makeExampleData(what="gs", gnames=names(se))
# (3) compiling artificial regulatory network 
grn <- makeExampleData(what="grn", nodes=names(se))
# (4) plot consistency graph
ggeaGraph(gs=gs[[1]], grn=grn, se=se)
# (5) get legend
ggeaGraphLegend()



# currently supported methods
nbeaMethods()

# (3) make 2 artificially enriched sets:
sig.genes <- names(se)[rowData(se)$ADJ.PVAL < 0.1]
gs[[1]] <- sample(sig.genes, length(gs[[1]])) 
gs[[2]] <- sample(sig.genes, length(gs[[2]]))   

# (4) gene regulatory network 
grn <- makeExampleData(what="grn", nodes=names(se))
    
# (5) performing the enrichment analysis
ea.res <- nbea(method="ggea", se=se, gs=gs, grn=grn)

# (6) result visualization and exploration
gsRanking(ea.res, signif.only=FALSE)

# using your own function as enrichment method
dummyNBEA <- function(se, gs, grn)
{
  sig.ps <- sample(seq(0,0.05, length=1000),5)
  insig.ps <- sample(seq(0.1,1, length=1000), length(gs)-5)
  ps <- sample(c(sig.ps, insig.ps), length(gs))
  score <- sample(1:100, length(gs), replace=TRUE)
  res.tbl <- cbind(score, ps)
  colnames(res.tbl) <- c("SCORE", "PVAL")
  rownames(res.tbl) <- names(gs)
  return(res.tbl[order(ps),])
}

ea.res2 <- nbea(method=dummyNBEA, se=se, gs=gs, grn=grn)
gsRanking(ea.res2)



hsa.grn <- compileGRN(org="hsa", db="kegg")
head(hsa.grn)

nbea.res <- nbea(method="ggea", se=allSE, gs=hsa.gs, grn=hsa.grn)
gsRanking(nbea.res)
