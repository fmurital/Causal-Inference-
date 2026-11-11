# Journal Paper

**Smoking During Pregnancy and Preterm Birth: A Doubly Robust Analysis of U.S. Births**

Faruk Muritala and Dhrubajyoti Ghosh, Ph.D., Kennesaw State University.

This folder holds the manuscript in preparation for journal submission, with the code, tables and figures it is built from. It builds on the SEASUG (CSOL) 2026 paper in [`../1_CSOL_SEASUG_2026_Presentation`](../1_CSOL_SEASUG_2026_Presentation) and the posters in [`../2_Conference_Posters`](../2_Conference_Posters).

| Path | Contents |
|---|---|
| `paper/manuscript.pdf` | The current draft (compiled) |
| `paper/manuscript.tex` | LaTeX source; numbers are inserted from `paper/gen/numbers.tex`, so no figure in the text is typed by hand |
| `paper/gen/` | Numbers and table fragments generated from the result tables by `paper/make_numbers.py` |
| `paper/figs/` | Figures used in the manuscript |
| `paper/preamble.tex`, `paper/refs_ms.tex` | Page setup and reference list |
| `code/` | `Smoking_Preterm_Full_Analysis.R`, the single script that produces every number and figure |
| `results/tables/` | Result tables from the full run (CSV) |
| `UPDATES.md` | What the draft adds to the CSOL paper, and what is still open |

## Rebuild the paper

1. Upload the contents of `paper` to a new Overleaf project and compile `manuscript.tex` with pdfLaTeX (the default). No extra packages are needed.
2. To regenerate the numbers after a new run, copy the result tables into `results/tables`, then run `python3 paper/make_numbers.py` from the `paper` folder.

## Status

Draft for review by the authors. The trimester-of-smoking section is marked in red in the manuscript and will be completed when the trimester data file is added to the `data` folder.
