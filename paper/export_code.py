"""Build the code appendix (code.tex): current sources in full, plus the diffs
that introduced and removed each defect the paper describes.

Full listings of all 17 recovered solver versions would run to thousands of
lines and obscure the point. Instead this gives every current source in full,
then for each issue a unified diff between the version before and the version
after, which is exactly the claim the paper makes about that issue. The
complete set of versions is in paper/variants/ for anyone who wants to run one.

Usage: python export_code.py
"""
import difflib
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, ".."))
VAR = os.path.join(HERE, "variants")

SOURCES = [
    ("qp.py", "The search: bounds, pruning and the two exact leaf solvers."),
    ("run_parallel.py", "Cuts the tree into independent subtrees and runs them."),
    ("split_task.py", "Splits one subtree further, for the case of Issue 7."),
    ("audit_results.py", "Audits a finished run for complete coverage."),
    ("test_qp.py", "Cross-check against brute force over odd squares."),
    ("test_leaves.py", "Planted solutions for the two closing lemmas."),
    ("peek.py", "Samples which states the search is visiting."),
    ("qfast.py", "Compiled search that closed the last omega = 17 node."),
    ("qfast_verify.py", "Settles qfast survivors with the exact search."),
    ("run_qfast.py", "Runs qfast over the unresolved node."),
    ("test_qfast.py", "Planted-solution test for qfast."),
    ("spot_qrange.py", "Spot check of the qfast node by the exact search."),
    ("qnode.py", "Compiled search of any node with two primes left (omega = 7)."),
    ("run_qnode.py", "Full omega = K run with the node kernel."),
    ("qchunks.py", "Splits dear nodes by the range of q."),
    ("finalize_omega7.py", "Combines the omega = 7 records and checks coverage."),
    ("test_qnode.py", "Node kernel against brute-force solutions and omega <= 6."),
]

# The conversation exporter is listed only in the private build; the public
# repository does not carry the conversation (QP_PUBLIC=1).
TOOLS = [
    ("paper/reconstruct_variants.py", "Rebuilds every source version from the transcript."),
    ("paper/verify_variants.py", "Runs each faulty version against its own defect."),
    ("paper/export_code.py", "Builds this appendix."),
]

if not os.environ.get("QP_PUBLIC"):
    TOOLS.insert(2, ("paper/export_conversation.py", "Builds the conversation appendix."))

DIFFS = [
    ("qp.v01.py", "qp.v07.py", "Issue 1: incomplete root finding, and its fix",
     "Before, the last-prime leaf bisects only the stretch beyond the maximum of "
     "$P_a$ and records an unresolved state for the rest. After, poly\\_roots "
     "locates the turning points from the single sign change of $P_a''$ and "
     "bisects every monotone stretch. Measured effect: 45 of 200 planted cases "
     "recovered and 190 flagged, against 200 recovered and none flagged."),
    ("qp.v14.py", "qp.v16.py", "Issue 8: the rule that does not terminate",
     "The added condition tests only \\texttt{hi}, which for an open component "
     "uses the limit $p/(p-1)$ and so does not change as the exponent bound "
     "rises. On the all-open path the condition therefore holds forever. "
     "Measured effect: an open exponent driven past $2 \\cdot 400$ after "
     "256,506 states."),
    ("qp.v16.py", "qp.v17.py", "Issue 8: the fix",
     "The condition now also tests the lower end of the admissible range, which "
     "does rise with the exponents, so it eventually fails. Measured effect: "
     "$10^6$ states on the same subtree with the deepest open exponent never "
     "exceeding $2 \\cdot 15$."),
]

HEAD = r"""\documentclass[9pt]{extarticle}
\usepackage[margin=0.75in]{geometry}
\usepackage{fontspec}
\setmonofont{DejaVuSansMono}[Extension=.ttf, Scale=0.82]
\usepackage{fancyvrb,fvextra,xcolor,hyperref}
\hypersetup{hidelinks}
\definecolor{add}{RGB}{20,110,40}
\definecolor{del}{RGB}{160,30,30}
\definecolor{hdr}{RGB}{90,90,90}
\newcommand{\difftitle}[1]{\par\medskip\noindent\textbf{#1}\par\smallskip}
\title{Appendix: source code}
\author{Ansel Marsh (an instance of Claude)}
\date{}
\begin{document}
\maketitle
\noindent Every source file as it currently stands, followed by the change that
introduced and the change that removed each defect discussed in the paper. The
seventeen recovered versions of the search itself are not listed in full --- they
would run to several thousand near-identical lines --- but all of them, with a
manifest of timestamps and hashes, are in \texttt{paper/variants/}, and each can
be run directly.
\tableofcontents
"""


def verbatim(text):
    # Verbatim needs no escaping; only guard against the end delimiter appearing.
    text = text.replace(r"\end{Verbatim}", r"\end {Verbatim}")
    return ("\\begin{Verbatim}[breaklines,breakanywhere,numbers=left,"
            "numbersep=4pt,fontsize=\\small]\n" + text.rstrip() + "\n\\end{Verbatim}\n")


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def diff_block(a, b):
    lines = list(difflib.unified_diff(
        read(os.path.join(VAR, a)).splitlines(),
        read(os.path.join(VAR, b)).splitlines(),
        fromfile=a, tofile=b, lineterm="", n=2))
    out = ["\\begin{Verbatim}[breaklines,breakanywhere,commandchars=\\\\\\{\\},"
           "fontsize=\\small]"]
    for L in lines:
        L = L.replace("\\", "\\textbackslash{}").replace("{", "\\{").replace("}", "\\}")
        if L.startswith("+++") or L.startswith("---") or L.startswith("@@"):
            out.append(f"\\textcolor{{hdr}}{{{L}}}")
        elif L.startswith("+"):
            out.append(f"\\textcolor{{add}}{{{L}}}")
        elif L.startswith("-"):
            out.append(f"\\textcolor{{del}}{{{L}}}")
        else:
            out.append(L)
    out.append("\\end{Verbatim}")
    return "\n".join(out) + "\n"


def main():
    parts = [HEAD, r"\section{Current sources}"]
    for name, blurb in SOURCES:
        parts += [f"\\subsection{{\\texttt{{{name.replace('_', chr(92) + '_')}}}}}",
                  blurb, verbatim(read(os.path.join(PROJ, name)))]
    parts.append(r"\section{Supporting tools}")
    for name, blurb in TOOLS:
        base = name.replace("_", chr(92) + "_")
        parts += [f"\\subsection{{\\texttt{{{base}}}}}", blurb,
                  verbatim(read(os.path.join(PROJ, name)))]
    parts.append(r"\section{The changes that introduced and removed each defect}")
    for a, b, title, blurb in DIFFS:
        parts += [f"\\subsection{{{title}}}", blurb, "", diff_block(a, b)]
    parts.append(r"\end{document}")
    with open(os.environ.get("QP_CODE_OUT", os.path.join(HERE, "code.tex")), "w", encoding="utf-8") as f:
        f.write("\n".join(parts) + "\n")
    print("wrote code.tex")


if __name__ == "__main__":
    main()
