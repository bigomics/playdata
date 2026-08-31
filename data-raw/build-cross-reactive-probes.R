## Build CROSS_REACTIVE_PROBES: Illumina methylation probes that map to more
## than one genome location, so their signal is not the CpG it is attributed
## to. Excluded from an EWAS by the "Cross-reactive" mask in the Methylome
## Profiler; see mp_xreactive_probes() in components/app_methylome.
##
## Why extdata and not data/: same reason as the EWAS catalog next door - a
## LazyData object stays in the namespace for the life of the session once
## touched, and this is read once per fit. As an .rds the reader owns it.
##
## Why here and not in the app repo: it is published reference data, not app
## code. It lived at components/app_methylome/inst/masking/ as a 2.7 MB CSV
## behind a .gitignore negation, and the loader tried a list of candidate
## paths because the app's working directory varies - which resolved to
## nothing once, silently, leaving 53,498 probes tested while the UI reported
## a mask had been applied. playdata::get_file() is system.file(mustWork =
## TRUE), so that failure mode is gone.
##
## SOURCE   Three published blacklists, carried in the `source` column:
##            chen2013_450k       29,233 rows  (450K)
##            mccartney2016_epic  44,210 rows  (EPIC)
##            zhou2017_mask       30,681 rows  (non-SNP general mask)
##          104,124 rows, 53,498 unique probes - a probe flagged by two lists
##          appears twice, and the column is kept so provenance survives.
##
## PROVENANCE GAP: the download URLs were never recorded. The CSV arrived
##          pre-assembled in omicsplayground 2825f3bc0 with only these three
##          labels, and no build script. The labels identify the papers, but
##          the exact files and accession dates are not known, so this script
##          converts what we have rather than rebuilding from source. Before
##          the next refresh, resolve each label to its published supplement
##          and rewrite this as a real download.
## RE-RUN   only after that gap is closed.

csv <- "../omicsplayground-methylome-app/components/app_methylome/inst/masking/cross_reactive_probes.csv"
stopifnot(file.exists(csv))

d <- utils::read.csv(csv, stringsAsFactors = FALSE)
stopifnot(identical(colnames(d), c("probe", "source")))

CROSS_REACTIVE_PROBES <- data.frame(
  probe  = factor(d$probe),
  source = factor(d$source),
  stringsAsFactors = FALSE
)
attr(CROSS_REACTIVE_PROBES, "sources") <- levels(CROSS_REACTIVE_PROBES$source)

message(nrow(CROSS_REACTIVE_PROBES), " rows, ",
        nlevels(CROSS_REACTIVE_PROBES$probe), " unique probes, from ",
        nlevels(CROSS_REACTIVE_PROBES$source), " lists")

saveRDS(CROSS_REACTIVE_PROBES, "inst/extdata/cross-reactive-probes.rds",
        compress = "xz")
