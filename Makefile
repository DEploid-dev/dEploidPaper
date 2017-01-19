.PHONY: all clean
all: clean natMethods.pdf natMethodsSupplement.pdf

mainfigures = $(shell grep png bioInfo.tex | sed -e "s/^.*{/figures\//g" -e "s/\}//g" )
supDEploidfigures = $(shell grep png bioInfoSupplementDEploid.tex | sed -e "s/^.*{//g" -e "s/\}//g" )
supCoveragefigures = $(shell grep png bioInfoSupplementCoverage.tex | sed -e "s/^.*{//g" -e "s/\}//g" )
suptex = $(shell grep "\.tex" bioInfoSupplement.tex | sed -e "s/^.*{//g" -e "s/\}//g" )

coverLetter.pdf: coverLetter.tex
	pdflatex coverLetter.tex

bioInfo.pdf: bioInfo.tex
	pdflatex bioInfo.tex
	pdflatex bioInfo.tex

natMethods.pdf: natMethods.tex
	pdflatex natMethods.tex
	pdflatex natMethods.tex

natMethodsSupplement.pdf: natMethodsSupplement.tex
	pdflatex natMethodsSupplement.tex
	pdflatex natMethodsSupplement.tex

natMethods.bbl: natMethods.bib
	pdflatex natMethods.tex
	bibtex natMethods.aux

bioInfoSupplement.pdf: bioInfoSupplement.tex ${supDEploidfigures} ${supCoveragefigures} ${suptex} supplementReset.tex
	pdflatex bioInfoSupplement.tex
	pdflatex bioInfoSupplement.tex

clean:
	rm -f *.blg *snm *nav *.bbl *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out bioInfoSupplement.pdf natMethods.pdf natMethodsSupplement.pdf

plain.pdf: bioInfo.tex Makefile
	sed -e "s/bioinfo/article/" \
	 -e "s/\\\usepackage{todonotes}/\\\usepackage{fontenc,inputenc,crop,graphicx,amsmath,array,color,amssymb,flushend,stfloats,amsthm,chngpage,times,fullpage, amsmath,natbib,todonotes}\n\\\def\\\address#1{\\\global\\\def\\\@issue{#1}}\n\\\def\\\history#1{\\\global\\\def\\\@history{#1}}\n\\\def\\\abstract#1{\\\global\\\def\\\@abstract{#1}}\n\\\def\\\corresp#1{\\\global\\\def\\\@corresp{#1}}\n/" \
	 -e "/copyrightyear/d" \
	 -e "s/\\\title\[Deconvoluting multiple infections\]/\\\title/" \
	 -e "s/\\\author\[Zhu \\\textit{et~al}.\]/\\\author/" \
	 -e "/Advance Access/d" \
	 -e "/Original Paper/d" \
	 -e "/\\\firstpage{1}/d" \
	 -e "/subtitle/d" \
	 -e "s/\\\sfb//g" \
	 -e "s/\\\sf//g" \
	 -e "/\\\history{/d" \
	 -e "/\\\maketitle/d" \
	 -e "s/\\\editor{Associate Editor: \\\textcolor{red}{XXXXXXX}}/\\\maketitle\n\\\noindent1. Wellcome Trust Centre for Human Genetics, University of Oxford, Oxford, UK \\\\\\\2. Medical Research Council (MRC) Centre for Genomics and Global Health, University of Oxford, Oxford, UK  \\\\\\\3. Wellcome Trust Sanger Institute, Hinxton, UK \\\\\\\4. Big Data Institute, Li Ka Shing Centre for Health Information and Discovery, University of Oxford, Oxford, UK \\\\\\\{*\} To whom correspondence should be addressed.\\\\\\\/" \
	 -e "s/\\\abstract{\\\textbf{Motivation:}/\\\begin{abstract}\\\\\\\\\\\noindent\\\textbf{Motivation:}  \\\\\\\ \\\noindent /g" \
	 -e "s/\\\textbf{Supplementary information:} Supplementary data are available at \\\textit{Bioinformatics} online.}/\\\end{abstract}/g" \
	 -e "s/\\\textbf{Results:}/\\\textbf{Results:}  \\\\\\\ \\\noindent /g" \
	 -e "s/\\\textbf{Availability and implementation:}/\\\textbf{Availability and implementation:}  \\\\\\\ \\\noindent /g" \
	 -e "s/\\\textbf{Contact:}/\\\textbf{Contact:}  \\\\\\\ \\\noindent /g" \
	 -e "s/0.45\\\textwidth/0.9\\\textwidth/g" \
	 -e "s/0.5\\\textwidth/0.9\\\textwidth/g" \
	 -e "s/.45\\\textwidth/0.9\\\textwidth/g" \
	  bioInfo.tex > plain.tex
	pdflatex plain.tex
	pdflatex plain.tex

pf3kFigures = $(shell grep "\.tex" pf3kDEploidNotes.tex | sed -e "s/^.*{//g" -e "s/\}//g" )
pf3kDEploidNotes.pdf: pf3kDEploidNotes.tex ${pf3kFigures}
	pdflatex pf3kDEploidNotes.tex
	pdflatex pf3kDEploidNotes.tex
