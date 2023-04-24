#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `CCLE_RNA_CANCERTYPE`
#' @source unknown
"CCLE_RNA_CANCERTYPE"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `CCLE_RNA_CELLINE`
#' @source unknown
"CCLE_RNA_CELLINE"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `DICE_SIGNATURE1000`
#' @source unknown
"DICE_SIGNATURE1000"

#' Example pgx object
#'
#' This is an example pgx object already created. It comes from
#' Geiger2016-arginine.
#'
#' @format ## `GEIGER_PGX`
#' A pgx object
#' \describe{
#'   \item{name}{PGX name}
#'   \item{description}{PGX description}
#'   \item{date}{PGX data created}
#'   ...
#' }
#' @source Geiger2016
"GEIGER_PGX"

#' Some kind of sparse matrix dataset?
#'
#' No idea what this data is but it is used in `test_genesets`.
#'
#' @format ## `GSET_SPARSEG_XL`
#' A dgCMatrix / Matrix object
#' @source unknown
"GSET_SPARSEG_XL"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `GTEX_RNA_TISSUE_TPM`
#' @source unknown
"GTEX_RNA_TISSUE_TPM"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `HPA_RNA_CELLINE`
#' @source unknown
"HPA_RNA_CELLINE"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `IMMPROT_SIGNATURE1000`
#' @source unknown
"IMMPROT_SIGNATURE1000"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `IMMUNOSTATES_MATRIX`
#' @source unknown
"IMMUNOSTATES_MATRIX"


#' A dataset used
#'
#' No idea what this data is
#'
#' @format ## `L1000_ACTIVITYS_N20D1011`
#' @source unknown
"L1000_ACTIVITYS_N20D1011"

#' A dataset used
#'
#' No idea what this data is
#'
#' @format ## `L1000_REPRURPOSING_DRUGS`
#' @source unknown
"L1000_REPRURPOSING_DRUGS"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `LM22`
#' @source unknown
"LM22"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `RNA_TISSUE_MATRIX`
#' @source unknown
"RNA_TISSUE_MATRIX"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `GSETS`
#' @source unknown
"GSETS"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `iGSETS`
#' @source unknown
"iGSETS"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `GSET_GENES`
#' @source unknown
"GSET_GENES"

#' Get path to omp example dataset(s)
#'
#' `playdata` comes bundled with a number of sample files in its `inst/extdata`
#' directory. This function make them easy to access. This function was
#' taken from tidyverse/readr.
#'
#' @param file string. Name of file. If `NULL`, the example files will
#'   be listed.

#' @examples
#' example_file()
#' example_file("counts.csv")
#' @export
example_file <- function(file = NULL) {
    if (is.null(file)) {
        dir(system.file("extdata", package = "playdata"))
    } else {
        system.file("extdata", file, package = "playdata", mustWork = TRUE)
    }
}
