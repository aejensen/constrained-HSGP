#!/usr/bin/env python3
"""Export the editing sources as independent submission TeX files.

LaTeX and BibTeX resolve references in a temporary directory. The exported
sources contain ordinary theorem/proof environments, embedded bibliographies,
and (for the supplement) embedded main-paper reference values. Only figures
remain external. --pdf also keeps the resulting PDFs. --check does not compile.
"""
import argparse
from collections import defaultdict
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

VARIANTS = {"main": "main", "supplement": "supplement",
            "end": "proofs_at_end", "inline": "proofs_inline"}
INPUT = re.compile(r"\\input\{([^}]+)\}")
PAIR = re.compile(
    r"\\begin\{theoremEnd\}(?:\[([^]]+)\])?\{([^}]+)\}\[([^]]+)\]"
    r"(.*?)\\end\{theoremEnd\}"
    r"(?:\s*\\begin\{proofEnd\}(.*?)\\end\{proofEnd\})?", re.S)
LABEL = re.compile(r"\\label\{([^}]+)\}")
BIBLIOGRAPHY = re.compile(r"\\bibliographystyle\{[^}]+\}\s*\\bibliography\{[^}]+\}")
FIGURE = re.compile(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}")
COUNTERS = ("equation", "figure", "table", "theorem", "corollary",
            "proposition", "lemma", "assumption", "definition", "algorithm")


def select_layout(text, inline):
    """Resolve the manuscript's layout switches before following file inputs."""
    tokens = re.compile(r"\\ifdefined\s*\\(CombinedInlineProofs|CombinedStandaloneSupplement)"
                        r"|\\else\b|\\fi\b")
    stack, output, cursor = [], [], 0
    for match in tokens.finditer(text):
        if all(stack):
            output.append(text[cursor:match.start()])
        if match.group(1):
            stack.append(inline if match.group(1) == "CombinedInlineProofs" else False)
        elif not stack:
            raise ValueError("Unmatched layout conditional: " + match.group())
        elif match.group() == r"\else":
            stack[-1] = not stack[-1]
        else:
            stack.pop()
        cursor = match.end()
    if stack:
        raise ValueError("Unclosed manuscript layout conditional")
    output.append(text[cursor:])
    return "".join(output)


def expand(source, name, inline=False, active=()):
    """Inline local inputs, retaining the editing files unchanged."""
    if name == "proof_configuration.tex":
        return r"\renewcommand*{\theHequation}{\arabic{equation}}" + "\n"
    path = (source / name).resolve()
    if not path.suffix:
        path = path.with_suffix(".tex")
    if path in active:
        raise ValueError("Circular TeX input: " + str(path))
    text = select_layout(path.read_text(), inline)
    return INPUT.sub(lambda m: expand(source, m[1], inline, (*active, path)), text)


def remove_checks(text):
    """The submission has no development logging commands."""
    text = re.sub(r"\\CombinedResultCheck\{[^}]+\}\{[^}]+\}", "", text)
    return re.sub(r"\\CombinedProofCheck\{[^}]+\}", "", text)


def convert_results(text, inline=False):
    """Replace proof-routing environments and retain principal statement/proof pairs."""
    deferred = defaultdict(list)

    def convert(match):
        options, kind, title, statement, proof = match.groups()
        statement = remove_checks(statement)
        proof = remove_checks(proof) if proof is not None else None
        principal = options is not None and "combined main" in options
        if principal:
            category = re.search(r"category=([^,]+)", options).group(1)
            label = LABEL.search(statement).group(1)
            deferred[category].append(dict(kind=kind, title=title, label=label,
                                           statement=statement, proof=proof))
        result = f"\\begin{{{kind}}}[{title}]{statement}\\end{{{kind}}}"
        if proof is not None and (inline or not principal):
            result += "\n\\begin{proof}" + proof + "\\end{proof}"
        return result

    converted = PAIR.sub(convert, text)
    if r"\begin{theoremEnd}" in converted or r"\begin{proofEnd}" in converted:
        raise ValueError("Unrecognized theorem/proof source structure")
    return converted, deferred


def groups(text):
    """Read adjacent braced fields of a LaTeX auxiliary label, including nested braces."""
    result, start, depth = [], None, 0
    for i, char in enumerate(text):
        slashes = 0
        j = i - 1
        while j >= 0 and text[j] == "\\":
            slashes += 1
            j -= 1
        if slashes % 2:
            continue
        if char == "{":
            if depth == 0:
                start = i + 1
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                result.append(text[start:i])
            if depth < 0:
                raise ValueError("Unbalanced auxiliary label")
    if depth:
        raise ValueError("Unbalanced auxiliary label")
    return result


def read_labels(aux):
    labels = {}
    for line in aux.splitlines():
        if line.startswith(r"\newlabel{"):
            key, value = groups(line[len(r"\newlabel"):])
            labels[key] = value
    return labels


def continuation(labels):
    """Continue numbered objects after the main paper, as in the combined layout."""
    maxima = dict.fromkeys(COUNTERS, 0)
    for key, value in labels.items():
        if not key.endswith("@cref"):
            continue
        match = re.match(r"\{\[([^]]+)\]\[(\d+)\]", value)
        if match and match[1] in maxima:
            maxima[match[1]] = max(maxima[match[1]], int(match[2]))
    return "\n".join(f"\\setcounter{{{kind}}}{{{number}}}" for kind, number in maxima.items())


def imported_labels(labels, local):
    """Embed reference values, without importing main-paper citations or dead PDF links."""
    lines = []
    for key, value in labels.items():
        if key.removesuffix("@cref") in local:
            continue
        if not key.endswith("@cref"):
            fields = groups(value)
            if len(fields) >= 4:
                fields[3] = ""  # The destination belongs to the other PDF.
            value = "".join("{" + field + "}" for field in fields)
        lines.append(r"\expandafter\def\csname r@" + key + r"\endcsname{" + value + "}")
    return "\n".join(lines)


def restatement(record, labels):
    """Restate the exact statement with its main-paper number and equation tags."""
    number = groups(labels[record["label"]])[0]

    def equation(match):
        body = match[1]
        keys = LABEL.findall(body)
        if len(keys) != 1:
            raise ValueError("A restated numbered equation must have one label")
        tag = groups(labels[keys[0]])[0]
        return r"\begin{equation*}\tag{" + tag + "}" + LABEL.sub("", body) + r"\end{equation*}"

    statement = re.sub(r"\\begin\{equation\}(.*?)\\end\{equation\}",
                       equation, record["statement"], flags=re.S)
    statement = LABEL.sub("", statement)
    return ("\\begingroup\n\\renewcommand{\\RestatedHeading}{"
            + record["kind"].capitalize() + " " + number + "}\n"
            + "\\begin{restatedresult}[" + record["title"] + "]"
            + statement + "\\end{restatedresult}\n"
            + "\\begin{proof}[Proof of \\cref{" + record["label"] + "}]"
            + record["proof"] + "\\end{proof}\n\\endgroup\n")


def prepare(source, inline=False):
    """Read the common preamble, main body, and supplementary body once per layout."""
    preamble = expand(source, "preamble.tex", inline)
    preamble = re.sub(r"^\\usepackage(?:\[[^]]*\])?\{(?:proof-at-the-end|xr-hyper)\}\s*\n",
                      "", preamble, flags=re.M)
    preamble += ("\n\\newcommand{\\RestatedHeading}{}\n"
                 "\\newtheorem*{restatedresult}{\\RestatedHeading}\n")
    main, deferred = convert_results(expand(source, "main_content.tex", inline), inline)
    supplement, _ = convert_results(expand(source, "supplement_content.tex", inline), True)
    if not inline and r"\begin{proof}" in main:
        raise ValueError("A main-paper proof is not assigned to a deferred category")
    categories = re.findall(r"\\printProofs\[([^]]+)\]", supplement)
    if not inline and (len(categories) != len(set(categories)) or set(categories) != set(deferred)):
        raise ValueError("Supplementary proof categories do not match the main results")
    return preamble, main, supplement, deferred


def run(command, folder):
    result = subprocess.run(command, cwd=folder, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(" ".join(command) + " failed:\n" + (result.stdout + result.stderr)[-6000:])


def compile_tex(folder, name, bibtex=False):
    """Resolve citations/references and return the aux text. All products stay temporary."""
    command = ["pdflatex", "-shell-escape", "-interaction=nonstopmode",
               "-halt-on-error", name + ".tex"]
    aux = folder / (name + ".aux")
    for iteration in range(8):
        previous = aux.read_text() if aux.exists() else None
        run(command, folder)
        if bibtex and iteration == 0:
            run(["bibtex", name], folder)
            continue
        text = aux.read_text()
        log = (folder / (name + ".log")).read_text(errors="replace")
        if text == previous and not re.search(
                r"(?:Citation|Reference).*undefined|There were undefined|multiply defined|"
                r"Label\(s\) may have changed|Rerun to get", log):
            return text
    raise RuntimeError("References did not resolve in " + name)


def finish(folder, name, text):
    """Embed the document's own bibliography and verify the resulting single TeX file."""
    (folder / (name + ".tex")).write_text(text)
    aux = compile_tex(folder, name, bibtex=True)
    bbl = (folder / (name + ".bbl")).read_text()
    # Each document's citations select its own bibliography. Main citations
    # are never imported into the supplement with its reference labels.
    cited = {key.strip() for item in re.findall(r"\\citation\{([^}]+)\}", aux)
             for key in item.split(",")}
    printed = set(re.findall(r"\\bibitem(?:\[[^]]*\])?\{([^}]+)\}", bbl))
    if cited != printed:
        raise RuntimeError("Bibliography mismatch in " + name)
    text = BIBLIOGRAPHY.sub(lambda _: bbl, text)
    if INPUT.search(text) or r"\externaldocument" in text or BIBLIOGRAPHY.search(text):
        raise RuntimeError("An external text dependency remains in " + name)
    (folder / (name + ".tex")).write_text(text)
    # Delete aux/bbl to check compilation without either document's old files.
    for suffix in (".aux", ".bbl"):
        (folder / (name + suffix)).unlink(missing_ok=True)
    aux = compile_tex(folder, name)
    return text, read_labels(aux)


def export(source, output, targets, pdf=False):
    """Export requested layouts, using a temporary main build for supplement references."""
    if output.resolve() == source.resolve():
        raise ValueError("Use a separate output directory to preserve the editing sources")
    preamble, main, supplement, deferred = prepare(source)
    output.mkdir(parents=True, exist_ok=True)
    completed = {}
    with tempfile.TemporaryDirectory(prefix="hsgp_submission_") as temporary:
        folder = Path(temporary)
        shutil.copy2(source / "references.bib", folder / "references.bib")
        for figure in set(FIGURE.findall(preamble + main + supplement)):
            destination = folder / figure
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / figure, destination)
        main_text = preamble + "\n\\begin{document}\n" + main + "\n\\end{document}\n"
        if any(target in targets for target in ("main", "supplement", "end")):
            print("Resolving main-paper references and bibliography...", flush=True)
            completed["main"], labels = finish(folder, "main", main_text)
        for variant in targets:
            if variant == "main":
                continue
            print("Preparing " + variant + "...", flush=True)
            if variant == "inline":
                head, body, tail, _ = prepare(source, inline=True)
                text = head + "\n\\begin{document}\n" + body + "\n\\clearpage\n\\appendix\n" + tail
            else:
                tail = re.sub(r"\\printProofs\[([^]]+)\]",
                    lambda m: "\n".join(restatement(r, labels) for r in deferred[m[1]]), supplement)
                if variant == "end":
                    text = preamble + "\n\\begin{document}\n" + main + "\n\\clearpage\n\\appendix\n" + tail
                else:
                    local = set(LABEL.findall(tail))
                    text = (preamble + "\n" + imported_labels(labels, local)
                        + "\n\\title{\\CombinedArticleTitle\\\\[1ex]{\\large Supplementary material}}\n"
                        + "\\hypersetup{pdftitle={Supplementary material: \\CombinedArticleTitle}}\n"
                        + "\\begin{document}\n\\maketitle\n" + continuation(labels)
                        + "\n\\appendix\n" + tail
                        + "\n\\bibliographystyle{plainnat}\n\\bibliography{references}\n")
            text += "\n\\end{document}\n"
            completed[variant], _ = finish(folder, VARIANTS[variant], text)
        for variant in targets:
            name = VARIANTS[variant]
            (output / (name + ".tex")).write_text(completed[variant])
            if pdf:
                shutil.copy2(folder / (name + ".pdf"), output / (name + ".pdf"))
            for figure in set(FIGURE.findall(completed[variant])):
                destination = output / figure
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(folder / figure, destination)
            print(output / (name + ".tex"))


def main(argv=None, source=None, default_output=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=source or Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path, default=default_output or Path("submission"))
    parser.add_argument("--target", choices=("submission", "all", *VARIANTS), default="submission")
    parser.add_argument("--pdf", action="store_true", help="Also keep PDFs in the export directory")
    parser.add_argument("--check", action="store_true", help="Check source expansion without compiling or exporting")
    args = parser.parse_args(argv)
    if args.check:
        for inline in (False, True):
            preamble, body, tail, deferred = prepare(args.source, inline)
            if INPUT.search(preamble + body + tail):
                raise ValueError("An unexpanded input remains")
            for figure in set(FIGURE.findall(preamble + body + tail)):
                if not (args.source / figure).is_file():
                    raise FileNotFoundError(args.source / figure)
        if not (args.source / "references.bib").is_file():
            raise FileNotFoundError(args.source / "references.bib")
        print("Source expansion checked. No compilation or exported files.")
        return
    targets = list(VARIANTS) if args.target == "all" else (
        ["main", "supplement"] if args.target == "submission" else [args.target])
    export(args.source.resolve(), args.output.resolve(), targets, args.pdf)


if __name__ == "__main__":
    main()
