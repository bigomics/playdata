
#' Get path to playdata file in inst/extdata
#'
#' `playdata` comes bundled with a number of non-RDA files in its `inst/extdata`
#' directory. This function make them easy to access. This function was
#' taken from tidyverse/readr.
#'
#' @param file string. Name of file. If `NULL`, the example files will
#'   be listed.

#' @examples
#' get_file()
#' get_file("kinase_substrates_kea.gmt")
#' @export
get_file <- function(file = NULL) {
    if (is.null(file)) {
        dir(system.file("extdata", package = "playdata"))
    } else {
        system.file("extdata", file, package = "playdata", mustWork = TRUE)
    }
}
