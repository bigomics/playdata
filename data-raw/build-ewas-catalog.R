## Build EWAS_CATALOG: which traits each CpG has already been associated with,
## and in how many study analyses - the "is this hit known or new?" annotation
## on an EWAS hit list.
##
## Why bundled: the app used to ask ewascatalog.org/api/ at runtime, one HTTP
## GET per CpG. That needs outbound network from the container, is slow enough
## to need a button, and pins us to whatever is live. The catalog publishes the
## whole thing as a flat file, so we ship a version we control.
##
## Why extdata and not data/: a LazyData object stays in the package namespace
## for the life of the R session once touched - 129 MB resident in every Shiny
## worker, for a table read once per hit list. As an .rds the reader owns it and
## the GC takes it back; see playbase.epigenetics::ewas_catalog_lookup().
##
## SOURCE   https://www.ewascatalog.org/download/
## VERSION  the results/studies pair Last-Modified 2026-04-09 (the catalog
##          publishes no version string; the mtime is the release marker)
## RE-RUN   when the catalog publishes a new release (roughly yearly)
##
## Reduction: results is 174 MB gzipped / 8.0M rows of per-association summary
## statistics. Nothing but the CpG and the study it came from is ever rendered,
## so beta/se/p/position/gene are dropped, the study id is resolved to its trait
## through the studies metadata file, and the pairs are counted. 8.0M rows of 13
## columns collapse to 6.3M rows of 3.

base <- "https://www.ewascatalog.org/static//docs/"
options(timeout = 3600) # the results file is 174 MB; the 60s default never finishes

fetch <- function(f) {
  p <- file.path(tempdir(), f)
  if (!file.exists(p)) {
    message("downloading ", f, " ...")
    utils::download.file(paste0(base, f), p, quiet = TRUE, mode = "wb")
  }
  p
}

## study_id -> trait. study_id is one analysis within a paper ("22232023_smoking
## _current_vs_never_smoking"), which is the unit the old API counted, so
## counting distinct study_ids reproduces the numbers the app used to show.
studies <- data.table::fread(fetch("ewascatalog-studies.txt.gz"),
                             select = c("study_id", "trait"), sep = "\t",
                             quote = "", colClasses = "character")

## Column 1 is cpg and column 6 is study_id; the file repeats the cpg column at
## position 7, so select by index rather than by name.
res <- data.table::fread(fetch("ewascatalog-results.txt.gz"),
                         select = c(1L, 6L), col.names = c("cpg", "study_id"),
                         sep = "\t", quote = "", colClasses = "character")

res <- unique(res) # one CpG can appear twice in an analysis
res[studies, trait := i.trait, on = "study_id"]
res <- res[!is.na(trait) & nzchar(trait)]

EWAS_CATALOG <- res[, .(n = data.table::uniqueN(study_id)), by = .(cpg, trait)]
data.table::setorder(EWAS_CATALOG, cpg, -n, trait)
EWAS_CATALOG <- data.frame(
  cpg = factor(EWAS_CATALOG$cpg),
  trait = factor(EWAS_CATALOG$trait),
  n = as.integer(EWAS_CATALOG$n),
  stringsAsFactors = FALSE
)
attr(EWAS_CATALOG, "source") <- "https://www.ewascatalog.org/download/"
attr(EWAS_CATALOG, "version") <- "2026-04-09"

message(nrow(EWAS_CATALOG), " CpG-trait pairs, ",
        nlevels(EWAS_CATALOG$cpg), " CpGs, ",
        nlevels(EWAS_CATALOG$trait), " traits")

saveRDS(EWAS_CATALOG, "inst/extdata/ewas-catalog.rds", compress = "xz")
