
<!-- README.md is generated from README.Rmd. Please edit that file -->

# playdata

<!-- badges: start -->
<!-- badges: end -->

The playdata package holds all of the internal datasets for the
`playbase` package.

# Adding Rda datasets

If you have a dataset to add, you should do the following:

- read the dataset into R memory via the console
- make sure the dataset variable is named as it will be called in
  playdata
- run `usethis::use_data(DATASET_NAME)`
- add a documentation entry in `R/playdata-data.R`
- run `devtools::document()` to update the documentation

Now you can call the dataset via `playdata::DATASET_NAME`

# Updating an Rda dataset

If you have a dataset to update, you should do the following:

- read the dataset into R memory via the console
- make sure the dataset variable is named as it will be called in
  playdata
- run `usethis::use_data(DATASET_NAME, overwrite = TRUE)`
- update documentation entry in `R/playdata-data.R`, if necessary
- run `devtools::document()` to update the documentation

# Adding non-Rda datasets

If you have a file to add that is not in Rda format, then do the
following:

- add it to the `inst/extdata` folder

Then you can get the path to the file via `playdata::get_file(FILENAME)`

# Accessing playdata datasets dynamically

If you want to access playdata dataset dynamically, where for example
the dataset from playdata to be accessed is only known via a variable
name, then use the function `playdata::get_data(...)`.

For example, if you want to access the dataset `BLUERED`, you can simply
call `playdata::BLUERED`. But if `BLUERED` was a variable - e.g.,
`mycolor <- 'BLUERED'` - then you can access the dataset via
`playdata::get_data(mycolor)`.
