#!/bin/sh
set -e
# https://posit-dev.github.io/air/
air format . inst/CITATION
Rscript -e 'devtools::document()'
pkg="ispso_$(sed '/Version/!d; s/.* //' DESCRIPTION).tar.gz"
R CMD build .
R CMD check $pkg
R CMD INSTALL $pkg
