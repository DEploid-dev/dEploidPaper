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
	sed -e "s/bioinfo/article/" \
	 -e "s/\\\usepackage{todonotes}/\\\usepackage{fontenc,inputenc,crop,graphicx,amsmath,array,color,amssymb,flushend,stfloats,amsthm,chngpage,times,fullpage, amsmath,natbib,todonotes}\n\\\def\\\address#1{\\\global\\\def\\\@issue{#1}}\n\\\def\\\history#1{\\\global\\\def\\\@history{#1}}\n\\\def\\\abstract#1{\\\global\\\def\\\@abstract{#1}}\n\\\def\\\corresp#1{\\\global\\\def\\\@corresp{#1}}\n/" \
	 -e "/copyrightyear/d" \
	 -e "s/\\\title\[Deconvolute mixed genomes\]/\\\title/" \
	 -e "s/\\\author\[Zhu \\\textit{et~al}.\]/\\\author/" \
	 -e "/Advance Access/d" \
	 -e "/Original Paper/d" \
	 -e "/\\\firstpage{1}/d" \
	 -e "/subtitle/d" \
	 -e "s/\\\sfb//g" \
	 -e "s/\\\sf//g" \
	 -e "/\\\history{/d" \
	 -e "/\\\editor{/d" \
	 -e "/\\\begin{methods}/d" \
	 -e "/\\\end{methods}/d" \
	 -e "s/\\\maketitle/\\\maketitle\n\\\noindent1. Wellcome Trust Centre for Human Genetics, University of Oxford, Oxford OX3 7BN, UK \\\\\\\2. Big data institute, Li Ka Shing Centre for Health Information and Discovery, University of Oxford, Oxford OX3 7BN, UK\\\\\\\\*. To whom correspondence should be addressed.\n/" \
	  bioInfo.tex > plain.tex
	pdflatex plain.tex
	pdflatex plain.tex
