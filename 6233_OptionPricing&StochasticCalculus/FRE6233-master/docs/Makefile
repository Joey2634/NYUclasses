SRC = ../x64/Debug
XML = $(wildcard $(SRC)/*.xml)
HTML = $(subst .xml,.html, $(notdir $(XML)))

%.html : $(SRC)/%.xml
	xsltproc --path $(SRC) $< -o $@

all: $(HTML)
	cp $(SRC)/*.css $(SRC)/index.html .
