# Stage 12 — manuscript and supplement edits

Scripts that insert the stage 10 and stage 11 results into the AHJ manuscript files with python-docx.

| script | edits |
|---|---|
| `edit_main.py` | Methods: two new subsections (calibration/clinical utility/precision; cell-type composition and pathway-level modeling). Results: two new subsections with four paragraphs. Discussion: one new paragraph and a rewritten Limitations sentence. Abstract Results: one added sentence. Data and code availability: sealed predictions. References 26-32 appended (no renumbering of 1-25). |
| `edit_supp.py` | Supplementary Methods: two new paragraphs. Table notes S5-S7. Figures S16-S18 embedded with legends before the Supplementary References. |

Both read the submitted .docx from a path set at the top of the script and write a `*_revised.docx`
alongside it; neither modifies the original. Figures are embedded as 600-dpi PNGs; the matching
600-dpi LZW TIFFs are produced by the stage 10 and 11 scripts for journal upload.
