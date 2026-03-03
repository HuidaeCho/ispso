#!/bin/sh
set -e
# curl -LsSf https://github.com/posit-dev/air/releases/latest/download/air-installer.sh | sh
air format . inst/CITATION
Rscript -e 'devtools::document()'
pkg="ispso_$(sed '/Version/!d; s/.* //' DESCRIPTION).tar.gz"
R CMD build .
R CMD check $pkg
R CMD INSTALL $pkg
