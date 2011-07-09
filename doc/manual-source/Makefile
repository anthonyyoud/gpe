#Copyright 2011 Anthony Youd/Newcastle University
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#------------------------------------------------------------------------------
OBJECT		= manual
BIBTEX		= true
LATEX		= latex
SOURCES		= $(OBJECT).tex introduction.tex quickstart.tex equations.tex \
						nondim.tex numerics.tex file_reference.tex

#---- IMPLICIT RULES ----------------------------------------------------------
.SUFFIXES: .tex .dvi .ps .pdf .ps.gz .bbl

.ps.ps.gz:
	gzip -c $*.ps > $*.ps.gz

.dvi.ps:
	dvips -o $*.ps $*.dvi

.ps.pdf:
	ps2pdf $*.ps
	cp $*.pdf ../

ifeq ($(LATEX), pdflatex)
# Using pdflatex to compile.
.tex.pdf:
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_pdf $*
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_pdf $*
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_pdf $*
	mv $*_pdf.pdf $*.pdf
	cp $*.pdf ../$*.pdf
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
	cp $*_dvi.bbl $*.bbl
else
# Using pdflatex to compile.
.tex.bbl:
	TEXINPUTS=./inputs:$(TEXINPUTS) $(LATEX) --jobname $*_pdf $*
	mv $*_pdf.pdf $*.pdf
	BSTINPUTS=./inputs:$(BSTINPUTS) bibtex $*_pdf
	cp $*_pdf.bbl $*.bbl
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
