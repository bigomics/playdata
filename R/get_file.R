
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


#' Get a playdata dataset dynamically from a variable name
#'
#' This function allows users to access playdata datasets when the dataset to
#' be accessed is only known dynamically or through a variable name
#'
#' @param name string. Name of dataset.
#' @examples
#' all.equal(get_dataset('BLUERED'), playdata::BLUERED)
#' @export
get_data <- function(name) {
    return(
        eval(rlang::parse_expr(glue::glue("playdata::{name}")))
    )
}
