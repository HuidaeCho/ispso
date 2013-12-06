# $Id: Makefile,v 1.6 2011/12/18 16:01:58 geni Exp $
LATEX=pdflatex -src-specials
BIBTEX=bibtex
VERSION=$(shell cat VERSION)

all: ispso.pdf ispso-$(VERSION).zip

# for VIM-LaTeX
pdf: ispso.pdf

ispso.pdf: ispso.tex VERSION
	sed 's/^\(\\date{.*Script Version\)[^,]*/\1 '`cat VERSION`'/' ispso.tex > ispso.tmp
	mv ispso.tmp ispso.tex
	$(LATEX) ispso
	$(BIBTEX) ispso
	$(LATEX) ispso
	$(LATEX) ispso

ispso-$(VERSION).zip: ispso.R funcs.R test.R ispso.pdf COPYING
	mkdir ispso-$(VERSION)
	cp $^ ispso-$(VERSION)
	zip -r $@ ispso-$(VERSION)
	rm -rf ispso-$(VERSION)

clean:
	rm -f ispso.aux ispso.bbl ispso.blg ispso.log ispso.toc
