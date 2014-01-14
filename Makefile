bibtex_example.pdf: bibtex_example.tex bibtex_example.bib
	pdflatex bibtex_example.tex
	bibtex bibtex_example
	pdflatex bibtex_example.tex
	pdflatex bibtex_example.tex
