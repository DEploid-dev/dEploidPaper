.PHONY: all clean
all: clean bioInfo.pdf bioInfoSupplement.pdf

figures = $(shell grep png *.tex | sed -e "s/^.*{//g" -e "s/\}//g" )

bioInfo.pdf: bioInfo.tex bioInfoSupplement.pdf
	pdflatex bioInfo.tex
	pdflatex bioInfo.tex

bioInfoSupplement.pdf: bioInfoSupplement.tex
	pdflatex bioInfoSupplement.tex
	pdflatex bioInfoSupplement.tex

clean:
	rm -f *.blg *snm *nav *.bbl *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out bioInfo.pdf bioInfoSupplement.pdf


#grep png *.tex | sed -e "s/^.*{//g" -e "s/\}//g"
