build: doc
	R -e "devtools::build()"


doc:
	R -e "devtools::document()"
	R -e "devtools::build_vignettes()"

check: 
	R -e "devtools::check()"

install: 
	R -e "devtools::document()"
	R CMD INSTALL .

