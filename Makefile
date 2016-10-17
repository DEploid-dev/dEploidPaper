.PHONY: all clean
all: clean bioInfo.pdf

mainfigures = $(shell grep png bioInfo.tex | sed -e "s/^.*{/figures\//g" -e "s/\}//g" )
supfigures = $(shell grep png bioInfoSupplement*.tex | sed -e "s/^.*{/figures\//g" -e "s/\}//g" )
suptex = $(shell grep "\.tex" bioInfoSupplement.tex | sed -e "s/^.*{//g" -e "s/\}//g" )

bioInfo.pdf: bioInfo.tex ${mainfigures}
	pdflatex bioInfo.tex
	pdflatex bioInfo.tex

bioInfoSupplement.pdf: bioInfoSupplement.tex ${supfigures} ${suptex} supplementReset.tex
	pdflatex bioInfoSupplement.tex
	pdflatex bioInfoSupplement.tex

clean:
	rm -f *.blg *snm *nav *.bbl *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out bioInfo.pdf bioInfoSupplement.pdf

plain.pdf: bioInfo.tex
#	sed -e "s/bioinfo/article/" -e "s/todonotes/fullpage, amsmath,natbib,todonotes" bioInfo.tex > plain.tex
	pdflatex plain.tex
	pdflatex plain.tex
