main.pdf: main.tex
	pdflatex $<
	bibtex main 
	pdflatex $<
	pdflatex $<

clean:
	rm -f main.aux
	rm -f main.log
	rm -f main.blg
	rm -f main.bbl
	rm -f main-blx.bib
	rm -f main.run.xml

spotless: clean
	rm -f main.pdf

.PHONY: clean spotless
