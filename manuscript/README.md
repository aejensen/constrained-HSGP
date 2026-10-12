# Manuscript editing and submission

Keep the individual `.tex` files for editing. `main.tex` is the existing
combined editing entry point. The `theory/` files keep each result beside
its proof, and `figures/` contains the application plots.

From this directory, generate the two submission sources with:

```sh
python3 export.py
```

This writes `submission/main.tex`, `submission/supplement.tex`, and the
figures they use. The main source contains no proofs. The supplement
restates the principal results and includes their proofs, supporting theory,
numerical methods, Stan material, and application details. Each file has its
own embedded bibliography containing only the references it cites.

The supplement embeds the main paper's reference values and continues its
numbering, so it compiles without a main-paper `.aux` file. Both exported
sources use ordinary theorem and proof environments and have no `\input`
dependencies. Keep their `figures/` directory beside them for submission.

The exporter needs Python 3.9 or newer, `pdflatex`, BibTeX, and the LaTeX
packages already used by the manuscript. It compiles in a temporary directory
to resolve reference numbers and format each bibliography, then checks the
exported source without its intermediate files. It does not change the
editing sources or the existing `main.pdf`. `--pdf` also copies the resulting
PDFs into the export directory. Compilation intermediates are removed.

```sh
python3 export.py --pdf                 # Submission pair with PDFs
python3 export.py --target all --pdf    # All four versions
python3 export.py --target inline      # Only the version with proofs inline
python3 export.py --check              # Source checks without compilation
```

`--target` accepts `main`, `supplement`, `end`, and `inline`, in addition to
`submission` and `all`. `--output DIRECTORY` selects the destination. Files
with the chosen output names are replaced. The combined files are named
`proofs_at_end.tex` and `proofs_inline.tex`. Regenerate the submission pair
after edits so the supplement's references match the current main paper.

Once exported, a file can be compiled on its own, for example:

```sh
latexmk -pdf -shell-escape supplement.tex
latexmk -c supplement.tex
```

BibTeX and the editing files are unnecessary for this final compilation.

In the development repository, use `python3 scripts/export_manuscript.py`
from the repository root instead. It stages the current authoritative
manuscript and theory files and writes to `output/submission/`, using the
same exporter and options. The public directory's exporter works independently
of that development repository.
