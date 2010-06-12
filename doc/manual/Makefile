# $Id$
#------------------------------------------------------------------------------

# Generic Makefile for compiling LaTeX documents.
#
# * OBJECT: Change this to the name of the main LaTeX file without extension.
# * BIBTEX: Set this to true if you are using BibTeX for your citations; set it
#           to false otherwise.
# * LATEX: Set this to latex or pdflatex depending on what output format you
#          want and what figure formats are included in the document.
# * SOURCES: These are the individual tex files making up the document.
#
# Usage:
# make 		- Compile the document.  Output depends on the LATEX and BIBTEX
#           variables.  Default is PS output using BibTeX.
# make twoup 	- Compile as PS but giving output in 2-up format.
# make book 	- Compile as PS but giving output suitable for printing in
# 		          booklet form.  Send to the printer as a duplex job.
# make tidy	- Remove any auxiliary files created during the make process.
# make clean 	- As for make tidy but in addition remove the final output
# 		          files.
#
#------------------------------------------------------------------------------
OBJECT		= manual

# Set this variable to true if you are using BibTeX for your citations; set it
# to false otherwise.
BIBTEX		= true

# Set this variable to latex if you are including only PS and EPS figures in
# your document; the resulting output will be a PS file.  Set it to pdflatex if
# you are only including PDF, PNG, or JPG images; the resulting output will be
# a PDF file.
LATEX		= latex

# Add list of files here to allow make to only build what is necessary.  Your
# master document should be listed first.
SOURCES		= $(OBJECT).tex introduction.tex equations.tex quickstart.tex \
						nondim.tex numerics.tex

#---- IMPLICIT RULES ----------------------------------------------------------
.SUFFIXES: .tex .dvi .ps .pdf .ps.gz .bbl

.ps.ps.gz:
	gzip -c $*.ps > $*.ps.gz

.dvi.ps:
	dvips -o $*.ps $*.dvi

.ps.pdf:
	ps2pdf $*.ps

ifeq ($(LATEX), pdflatex)
# Using pdflatex to compile.
.tex.pdf:
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_pdf $*
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_pdf $*
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_pdf $*
	mv $*_pdf.pdf $*.pdf
endif

.tex.dvi:
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_dvi $*
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_dvi $*
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_dvi $*
	mv $*_dvi.dvi $*.dvi

ifeq ($(LATEX), latex)
# Using latex to compile.
.tex.bbl:
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_dvi $*
	mv $*_dvi.dvi $*.dvi
	BSTINPUTS=./inputs:$(BSTINPUTS) bibtex $*_dvi
	cp  $*_dvi.bbl  $*.bbl
else
# Using pdflatex to compile.
.tex.bbl:
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_pdf $*
	mv $*_pdf.pdf $*.pdf
	BSTINPUTS=./inputs:$(BSTINPUTS) bibtex $*_pdf
	cp  $*_pdf.bbl  $*.bbl
endif

#------------------------------------------------------------------------------
ifeq ($(LATEX), latex)
# Using latex to compile.
default: ps
else
# Using pdflatex to compile.
default: pdf
endif

all: $(OBJECT).dvi ps

#------------------------------------------------------------------------------
ps: $(OBJECT).ps

pdf: $(OBJECT).pdf

twoup: $(OBJECT).ps 
	psnup -n2 -Pa4 $(OBJECT).ps $(OBJECT).2up.ps

book: ps
	psbook $(OBJECT).ps book.ps
	psnup -n2 -Pa4 book.ps $(OBJECT).book.ps
	rm -f book.ps

#------------------------------------------------------------------------------
tidy:
	rm -f *.aux *.toc *.log *.out  \
	      *.bbl *.blg *.idx *.lot *.lof \
	      $(OBJECT)_dvi.* $(OBJECT)_pdf.* */*.aux
clean: tidy
	rm -f $(OBJECT).dvi $(OBJECT).pdf $(OBJECT).ps

#---- DEPENDENCIES ------------------------------------------------------------

ifeq ($(BIBTEX), true)
  # With BibTeX
  ifeq ($(LATEX), latex)
    # Using latex to compile.
    $(OBJECT).dvi: $(OBJECT).bbl $(OBJECT).tex $(SOURCES)
  else
    # Using pdflatex to compile.
    $(OBJECT).pdf: $(OBJECT).bbl $(OBJECT).tex $(SOURCES)
  endif
else
  # Without BibTeX
  ifeq ($(LATEX), latex)
    # Using latex to compile.
    $(OBJECT).dvi: $(OBJECT).tex
  else
    # Using pdflatex to compile.
    $(OBJECT).pdf: $(OBJECT).tex
  endif
endif

$(OBJECT).bbl: $(OBJECT).tex references.bib $(SOURCES)
$(OBJECT).ps: $(OBJECT).dvi
$(OBJECT).ps.gz: $(OBJECT).ps
