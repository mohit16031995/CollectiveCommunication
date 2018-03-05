PAPER=paper
TEX = $(wildcard *.tex)
BIB = references.bib
FIGS = $(wildcard img/*.pdf img/*.eps)

.PHONY: all clean

$(PAPER).pdf: $(TEX) $(BIB) $(FIGS) 
	echo $(FIGS)
	pdflatex $(PAPER)
	bibtex $(PAPER)
	pdflatex $(PAPER)
	pdflatex $(PAPER)

clean:
	rm -f *.aux *.bbl *.blg *.log *.out $(PAPER).dvi $(PAPER).pdf

