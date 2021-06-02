
# Examples

PTXQC will provide you with a report in Html/PDF format, summarizing the quality of your data.

Please refer to the publication (currently under review; link is coming soon) for details on which quality metrics are available and how they are scored.

## Obtaining example data
If you want to generate a report yourself, example input data can be found on the **PRIDE** repository under the
following accession numbers: [PXD003133 (MBR)][1], [PXD003134 (Myco)][2], and [PXD000427 (TMT)][3].
This data corresponds to the MaxQuant data used in the publication.

The [PTXQC-Basic_Guide_for_R_users][4] vignette has R code snippets to automatically download the data. If you are not an R person and have Windows as operating system, just download the zip files from the PRIDE link above, and run PTXQC via drag'n'drop (as described in the [PTXQC-DragNDrop][5] vignette).

## Report Layout

The first page shows an overview heatmap, representing scores for each Raw file for a set of quality metrics.
The screenshots given here, a from a rather small run of just 5 files. You can easily analyse a study with more than 100 files, but the report will get bigger :)

![Overview Heatmap](./example_heatmap.png?raw=true "Overview heatmap showing quality criteria for each LC-MS file")
 
After that, there is a summary of MaxQuant parameters (i.e. which FASTA files were used, was Match-between-runs activated etc...).

In the remainder of the report you can follow up on each individual metric and explore the reason for failure. The Html output additionally contains a
helpful description of each metric in case you have not read the [publication][JPR_paper] or need a refresh on what the plot actually tells you.

For example, the alignment in MaxQuant via the Match-between-runs did not work optimal across all Raw files in this case,
since MaxQuant could not normalize the LC gradients at the start and end of the run:
![Alignment Performance](./example_MBRalignment.png?raw=true "Alignment of 5 raw files (the first file serves as reference here)")

Also, you could check how a new LC gradient influences the identifications over time:
![Identifications over Retention Time](./example_IDoverRT.png?raw=true "Identifications over Retention Time")

Find a full report as [html][example_html] or [pdf][example_pdf].

See the package vignettes for documentation on how to create and customize a report.

  [1]: http://www.ebi.ac.uk/pride/archive/projects/PXD003133
  [2]: http://www.ebi.ac.uk/pride/archive/projects/PXD003134
  [3]: http://www.ebi.ac.uk/pride/archive/projects/PXD000427
  [4]: https://github.com/cbielow/PTXQC/blob/master/vignettes/PTXQC-Basic_Guide_for_R_users.Rmd
  [5]: https://github.com/cbielow/PTXQC/blob/master/vignettes/PTXQC-DragNDrop.Rmd
  [JPR_paper]: https://doi.org/10.1021/acs.jproteome.5b00780
  [example_html]: http://htmlpreview.github.io/?https://github.com/cbielow/PTXQC/blob/master/inst/examples/report_v1.0.10__txt_5files_withMatch-100min.html
  [example_pdf]: report_v1.0.10__txt_5files_withMatch-100min.pdf