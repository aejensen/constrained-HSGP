.PHONY: all submission main supplement end inline clean

all:
	python3 export.py --target all --pdf

submission:
	python3 export.py --pdf

main supplement end inline:
	python3 export.py --target $@ --pdf

clean:
	latexmk -c main.tex
	rm -f main-pratend*.tex
	rm -rf _minted _minted-main
