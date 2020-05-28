check:
	Rscript -e "devtools::check()"

bioccheck:
	Rscript -e "BiocCheck::BiocCheck('.')"

install:
	Rscript -e "devtools::install(build_vignettes = TRUE, upgrade_dependencies = FALSE)"
