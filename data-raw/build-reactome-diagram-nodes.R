## Build inst/extdata/reactome_diagram_nodes.tsv.gz: where each molecule sits on
## every bundled Reactome diagram, and the ids it can be matched by.
##
## Why: the bundled SVGs (build-reactome-diagrams.R) draw their labels as
## outlines, so a renderer cannot find a node by its text. Reactome publishes the
## layout of each diagram separately, at the same coordinates as the SVG; with it
## a renderer recolours the node's own shape instead.
##
## One row per Chemical or Protein node: diagram, x, y, width, height, keys.
## `keys` is `;`-joined -- gene symbols and the UniProt accession for a protein,
## `ChEBI_<id>` for a chemical. Reactome draws chemicals in their charged form
## (L-argininium(1+)) while datasets record the neutral compound, so each ChEBI id
## is expanded to every conjugate acid/base and tautomer of it.
## ponytail: Complex/EntitySet/Gene nodes are skipped -- one colour for several
## molecules says nothing; add them if a member-summary rule is ever agreed.
##
## Run after build-reactome-diagrams.R (it reads the diagram ids bundled there).
## Fetches two small JSON files per diagram from Reactome's download area.

base <- "https://reactome.org/download/current/diagram"
extdata <- file.path("inst", "extdata")
D <- sub("\\.svg\\.gz$", "", list.files(file.path(extdata, "reactome-svg"),
                                        pattern = "\\.svg\\.gz$"))

## ChEBI charge-state classes: conjugate acid (6), conjugate base (7), tautomer (11)
rel_file <- tempfile(fileext = ".tsv.gz")
utils::download.file(
  "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/relation.tsv.gz",
  rel_file, quiet = TRUE
)
rel <- utils::read.delim(rel_file, colClasses = "character")
rel <- rel[rel$relation_type_id %in% c("6", "7", "11"), c("init_id", "final_id")]
comp <- igraph::components(igraph::graph_from_data_frame(rel, directed = FALSE))$membership
members <- split(names(comp), comp)
chebi_class <- function(id) {
  if (id %in% names(comp)) members[[as.character(comp[[id]])]] else id
}

fetch_json <- function(url) {
  for (attempt in 1:3) {
    out <- tryCatch(jsonlite::fromJSON(url, simplifyVector = FALSE),
                    error = function(e) NULL)
    if (!is.null(out)) return(out)
    Sys.sleep(attempt * 2)
  }
  NULL
}

rows <- list()
for (i in seq_along(D)) {
  d <- D[[i]]
  layout <- fetch_json(sprintf("%s/%s.json", base, d))
  graph <- fetch_json(sprintf("%s/%s.graph.json", base, d))
  if (is.null(layout) || is.null(graph)) {
    message("skipped ", d, " (layout unavailable)")
    next
  }
  entities <- stats::setNames(graph$nodes, vapply(graph$nodes, function(e) {
    as.character(e$dbId)
  }, character(1)))
  for (n in layout$nodes) {
    if (!n$renderableClass %in% c("Chemical", "Protein")) next
    e <- entities[[as.character(n$reactomeId)]]
    if (is.null(e) || is.null(e$identifier)) next
    keys <- if (identical(n$renderableClass, "Chemical")) {
      paste0("ChEBI_", chebi_class(e$identifier))
    } else {
      c(unlist(e$geneNames), e$identifier)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      diagram = d, x = n$prop$x, y = n$prop$y,
      width = n$prop$width, height = n$prop$height,
      keys = paste(unique(keys), collapse = ";")
    )
  }
  if (i %% 100 == 0) message(i, "/", length(D), " diagrams")
}

nodes <- do.call(rbind, rows)
con <- gzfile(file.path(extdata, "reactome_diagram_nodes.tsv.gz"), "w")
utils::write.table(nodes, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
message("wrote ", nrow(nodes), " nodes over ", length(unique(nodes$diagram)),
        " of ", length(D), " diagrams")
