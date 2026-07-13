## Build the Reactome native-diagram assets shipped in inst/extdata:
##   - reactome-svg/R-HSA-<id>.svg.gz : Reactome's own pre-rendered diagram SVGs
##   - reactome_diagram_map.tsv       : every human pathway -> the diagram id
##                                      that depicts it (self, or nearest ancestor
##                                      that owns a diagram)
##
## Why: reactome.org's live diagram exporter
## (/ContentService/exporter/diagram/<id>.svg) is behind a Cloudflare bot
## challenge and returns 403 to server-side requests. The static download area is
## NOT challenged, and it publishes the same pre-rendered SVGs in bulk. We mirror
## those and precompute the sub-pathway -> diagram resolution offline, so the app
## never calls reactome.org at runtime.
##
## Re-run when Reactome ships a new release (~quarterly).

base <- "https://reactome.org/download/current"
extdata <- file.path("inst", "extdata")
svgdir <- file.path(extdata, "reactome-svg")
dir.create(svgdir, showWarnings = FALSE, recursive = TRUE)

tmp <- tempfile(fileext = ".tgz")
message("downloading diagrams.svg.tgz (~260 MB) ...")
utils::download.file(file.path(base, "diagrams.svg.tgz"), tmp, quiet = TRUE)

## extract only the human diagrams (R-HSA-*.svg), gzip each into inst/extdata
raw <- tempfile("reactome-svg-")
dir.create(raw)
utils::untar(tmp, files = NULL, exdir = raw) # entries are flat: R-HSA-<id>.svg
svgs <- list.files(raw, pattern = "^R-HSA-.*\\.svg$", full.names = TRUE, recursive = TRUE)
message("gzipping ", length(svgs), " human diagram SVGs ...")
unlink(list.files(svgdir, full.names = TRUE))
for (f in svgs) {
  R.utils::gzip(f, destname = file.path(svgdir, paste0(basename(f), ".gz")),
                overwrite = TRUE, remove = FALSE)
}

## D = diagram ids we actually have an SVG for (the map must only target these)
D <- sub("\\.svg\\.gz$", "", list.files(svgdir, pattern = "\\.svg\\.gz$"))

## parent<TAB>child relations, human only -> child -> parents lookup
rel <- read.delim(file.path(base, "ReactomePathwaysRelation.txt"),
                  header = FALSE, colClasses = "character",
                  col.names = c("parent", "child"))
rel <- rel[startsWith(rel$parent, "R-HSA-") & startsWith(rel$child, "R-HSA-"), ]
parents <- split(rel$parent, rel$child)
nodes <- union(D, union(rel$parent, rel$child))

Dset <- new.env(hash = TRUE)
for (d in D) assign(d, TRUE, Dset)

## nearest ancestor (incl. self) that owns a diagram, via BFS up the hierarchy
resolve <- function(x) {
  seen <- new.env(hash = TRUE)
  queue <- x
  while (length(queue)) {
    n <- queue[1]
    queue <- queue[-1]
    if (!is.null(Dset[[n]])) return(n)
    if (!is.null(seen[[n]])) next
    assign(n, TRUE, seen)
    queue <- c(queue, parents[[n]])
  }
  NA_character_
}
diagram <- vapply(nodes, resolve, character(1))
map <- data.frame(id = nodes, diagram = diagram, stringsAsFactors = FALSE)
map <- map[!is.na(map$diagram), ]
map <- map[order(map$id), ]
write.table(map, file.path(extdata, "reactome_diagram_map.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message("wrote ", length(D), " diagrams and ", nrow(map), " map entries (",
        round(100 * nrow(map) / length(nodes)), "% coverage)")
