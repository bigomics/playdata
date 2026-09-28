## Build the two EPIC v2 reference tables playbase.epigenetics reads:
##
##   inst/extdata/epicv2-replicates.rds     which replicate probe to keep per CpG
##   inst/extdata/cross-reactive-probes.rds + the EPIC v2 cross-hybridising
##                                            list, source "peters2024_epicv2"
##
## Run from the package root, after build-cross-reactive-probes.R when that one
## is re-run (this script replaces only its own source's rows, so either order
## of re-runs converges).
##
## ---- epicv2-replicates.rds -------------------------------------------------
##
## WHY      EPIC v2 names every probe <id>_<design suffix> (cg00000029_TC21) and
##          measures ~5,200 ids with 2-10 replicate probes. Clocks, cell
##          references, the EWAS Catalog and the masks are all keyed by the bare
##          id, so playbase.epigenetics::collapse_epicv2_replicates() keeps one
##          replicate per id. Never the mean: replicates differ in design type
##          and in how well they track true methylation, so their average is a
##          probe nobody characterised.
##
## SHAPE    one row per probe whose bare id is replicated in the manifest
##          (11,616 probes over 5,222 ids, ordered by id then rank):
##            probe_id     full IlmnID            (factor)
##            cpg          bare id                (factor)
##            recommended  keep this one          (logical, exactly one per cpg)
##            evidence     Peters 2024 verdict    (factor, NA when none)
##          attr "manifest" = "20a1": the Illumina manifest release, as the
##          IlluminaHumanMethylationEPICv2anno.20a1.hg38 package names it.
##          playbase.epigenetics refuses a table whose manifest differs from its
##          annotation package, so a manifest bump needs a paired release here.
##          attr "source" names both inputs.
##
## WHICH    Illumina's manifest does NOT designate a preferred replicate: its
##          only replicate field is Rep_Num, a synthesis counter. The one
##          published per-probe evaluation is Peters et al. 2024 (below), which
##          compared each replicate against matched EPIC v1 and WGBS and labels
##          every replicated probe in "Rep_results_by_NAME". Ranked, first wins:
##            1 "Superior probe" / "Superior by WGBS"   best on both counts
##            2 "Best sensitivity"                      tracks methylation change
##            3 "Best precision"
##            4 no verdict, "Insufficient evidence", or a verdict that only the
##              probe-set MEAN is best ("... group mean ...") - averaging is
##              ruled out above, so these fall through to the tie-break
##            5 "Inferior probe" / "Inferior by WGBS"
##          Ties go to a probe the annotation package can place (two sets mix
##          a chr0 replicate with a mapped one), then the lowest replicate
##          number, then the probe id. This is DMRcate's filter.strategy =
##          "sensitivity" (same author), made deterministic: DMRcate breaks
##          ties with sample().
##
## ---- cross-reactive-probes.rds, source peters2024_epicv2 --------------------
##
## Peters 2024 BLAT-aligned every EPIC v2 probe; "CH_BLAT" = Y marks the 30,627
## in-silico cross-hybridising probes (47+ matching bases off target), and the
## paper adds the 6,889 probes Illumina could not place ("CHR" = chr0, 97%
## multimappers) to reach its "total of 37,346". Both are taken, matching the
## in-silico definition of the three lists already in the file. Ids are bare,
## and a replicated id is flagged by the replicate that survives the collapse
## above - the probe whose betas a dataset actually holds. SNP-affected probes
## are not added here: the Methylome app masks those from the manifest's own
## SNP columns per array.
##
## SOURCE   Peters TJ et al. (2024) Characterisation and reproducibility of the
##          HumanMethylationEPIC v2.0 BeadChip for DNA methylation profiling.
##          BMC Genomics 25:251. doi:10.1186/s12864-024-10027-5
##          Additional file 4 (augmented EPIC v2 manifest, 1.2 GB CSV).
##          Bioconductor's EPICv2manifest package carries the same table, but
##          only from Bioconductor 3.19 and via AnnotationHub; the published
##          supplement is the citable artefact, pinned by checksum below.
##          Probe universe: the installed IlluminaHumanMethylationEPICv2manifest;
##          placement: IlluminaHumanMethylationEPICv2anno.20a1.hg38 (both
##          Bioconductor 3.19, both install on R >= 4.2 / minfi 1.48).
## RE-RUN   when Illumina ships a new EPIC v2 manifest (new annotation package
##          name) or a revised replicate evaluation is published.

suppressPackageStartupMessages({
  stopifnot(requireNamespace("data.table"), requireNamespace("minfi"))
})
ANNO <- "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
stopifnot(requireNamespace(ANNO))
MANIFEST <- strsplit(ANNO, ".", fixed = TRUE)[[1]][2] # "20a1"
SUFFIX <- "_[TB][CO][12][0-9]+$"

PETERS_URL <- paste0(
  "https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-024-10027-5/",
  "MediaObjects/12864_2024_10027_MOESM4_ESM.csv"
)
PETERS_MD5 <- "259eb024cf0ab064eb348c376c40dbd2"

options(timeout = 3600) # 1.2 GB
f <- file.path(tempdir(), "peters2024_epicv2_manifest.csv")
if (!file.exists(f)) {
  message("downloading Peters 2024 Additional file 4 ...")
  utils::download.file(PETERS_URL, f, quiet = TRUE, mode = "wb")
}
stopifnot(unname(tools::md5sum(f)) == PETERS_MD5)
peters <- data.table::fread(
  f, select = c("IlmnID", "CHR", "CH_BLAT", "Rep_results_by_NAME")
)
peters <- as.data.frame(peters)
rownames(peters) <- peters$IlmnID

## ---- replicate table --------------------------------------------------------

## The universe is every probe minfi can return a beta for: the manifest
## package's type I and II probes (not the rs genotyping probes).
MANIFEST_PKG <- "IlluminaHumanMethylationEPICv2manifest"
stopifnot(requireNamespace(MANIFEST_PKG))
mf <- getExportedValue(MANIFEST_PKG, MANIFEST_PKG)
ids <- c(minfi::getProbeInfo(mf, type = "I")$Name, minfi::getProbeInfo(mf, type = "II")$Name)
stopifnot(all(grepl(SUFFIX, ids)), !anyDuplicated(ids))
bare <- sub(SUFFIX, "", ids)
is_rep <- bare %in% bare[duplicated(bare)]

## Probes the annotation package can place. It leaves out Illumina's chr0
## (unmapped) and chrM probes, so their betas never reach an annotated matrix.
library(ANNO, character.only = TRUE) # getAnnotation() reads the search path
placed <- ids %in% rownames(minfi::getAnnotation(get(ANNO)))

## Every replicated probe is one Peters evaluated.
stopifnot(all(ids[is_rep] %in% peters$IlmnID))

evidence <- peters[ids[is_rep], "Rep_results_by_NAME"]
evidence[!nzchar(evidence)] <- NA_character_
tier <- ifelse(
  grepl("^Superior (probe|by WGBS)$", evidence), 1L,
  ifelse(evidence %in% "Best sensitivity", 2L,
    ifelse(evidence %in% "Best precision", 3L,
      ifelse(grepl("^Inferior", evidence), 5L, 4L)
    )
  )
)

reps <- data.frame(
  probe_id = ids[is_rep], cpg = bare[is_rep], evidence = evidence,
  tier = tier, unplaced = !placed[is_rep],
  rep_num = as.integer(sub("^.*_[TB][CO][12]", "", ids[is_rep])),
  stringsAsFactors = FALSE
)
reps <- reps[order(reps$cpg, reps$tier, reps$unplaced, reps$rep_num, reps$probe_id), ]
reps$recommended <- !duplicated(reps$cpg)

stopifnot(
  all(table(reps$cpg[reps$recommended]) == 1),
  setequal(reps$cpg[reps$recommended], unique(reps$cpg))
)

EPICV2_REPLICATES <- data.frame(
  probe_id = factor(reps$probe_id),
  cpg = factor(reps$cpg),
  recommended = reps$recommended,
  evidence = factor(reps$evidence),
  row.names = NULL
)
attr(EPICV2_REPLICATES, "manifest") <- MANIFEST
attr(EPICV2_REPLICATES, "source") <- paste0(
  MANIFEST_PKG, " ", utils::packageVersion(MANIFEST_PKG), " and ",
  ANNO, " ", utils::packageVersion(ANNO), " (Illumina EPIC v2.0 manifest ",
  MANIFEST, "); replicate choice from Peters et al. 2024, BMC Genomics 25:251, ",
  "Additional file 4 (md5 ", PETERS_MD5, "), Rep_results_by_NAME; built ",
  Sys.Date()
)

message(
  nrow(EPICV2_REPLICATES), " replicate probes over ",
  nlevels(EPICV2_REPLICATES$cpg), " ids; recommended by verdict:"
)
print(table(reps$evidence[reps$recommended], useNA = "ifany"))

saveRDS(EPICV2_REPLICATES, "inst/extdata/epicv2-replicates.rds", compress = "xz")

## ---- cross-reactive list ----------------------------------------------------

survives <- !is_rep
survives[match(reps$probe_id, ids)] <- reps$recommended
flagged <- peters$IlmnID[peters$CH_BLAT %in% "Y" | peters$CHR %in% "chr0"]
xr <- sort(unique(bare[survives & ids %in% flagged]))
message(length(flagged), " flagged probes -> ", length(xr), " bare ids after the collapse")

old <- readRDS("inst/extdata/cross-reactive-probes.rds")
old <- old[old$source != "peters2024_epicv2", , drop = FALSE]
CROSS_REACTIVE_PROBES <- data.frame(
  probe = factor(c(as.character(old$probe), xr)),
  source = factor(c(as.character(old$source), rep("peters2024_epicv2", length(xr)))),
  stringsAsFactors = FALSE
)
attr(CROSS_REACTIVE_PROBES, "sources") <- levels(CROSS_REACTIVE_PROBES$source)

message(nrow(CROSS_REACTIVE_PROBES), " rows, ",
        nlevels(CROSS_REACTIVE_PROBES$probe), " unique probes, from ",
        nlevels(CROSS_REACTIVE_PROBES$source), " lists")

saveRDS(CROSS_REACTIVE_PROBES, "inst/extdata/cross-reactive-probes.rds",
        compress = "xz")
