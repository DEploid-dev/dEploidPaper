.PHONY: all clean
all: clean bioInfo.pdf bioInfoSupplement.pdf

mainfigures = $(shell grep png bioInfo.tex | sed -e "s/^.*{/figures\//g" -e "s/\}//g" )
supfigures = $(shell grep png bioInfoSupplement*.tex | sed -e "s/^.*{/figures\//g" -e "s/\}//g" )

bioInfo.pdf: bioInfo.tex bioInfoSupplement.pdf ${mainfigures}
	pdflatex bioInfo.tex
	pdflatex bioInfo.tex

bioInfoSupplement.pdf: bioInfoSupplement*.tex ${supfigures}
	pdflatex bioInfoSupplement.tex
	pdflatex bioInfoSupplement.tex

clean:
	rm -f *.blg *snm *nav *.bbl *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out bioInfo.pdf bioInfoSupplement.pdf


#grep png *.tex | sed -e "s/^.*{//g" -e "s/\}//g"
