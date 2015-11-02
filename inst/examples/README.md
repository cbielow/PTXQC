
# Examples

PTXQC will provide you with a report in PDF/HTML format, summarizing the quality of your data.

Please refer to the publication (currently under review; link is coming soon) for details on which quality metrics are available and how they are scored.

## Obtaining example data
If you want to generate a report yourself, example input data can be found on the PRIDE repository under the following accession numbers: [PXD003133] [1]. This data corresponds to the MaxQuant data used in the publication.
The 'PTXQC-Basic_Guide_for_R_users' vignette has R code snippets to automate this. If you are not an R person and have Windows as operating system, just download the zip files from the PRIDE link above, and run PTXQC via drag'n'drop (as described in the 'PTXQC-DragNDrop' vignette).

## Report Layout

The first page in the PDF will be a summary of MaxQuant parameters (i.e. which FASTA files were used, was Match-between-runs activated etc...).

The second page is an overview heatmap.
The screenshots given here, a from a rather small run of just 5 files. You can easily analyse a study with more than 100 files, but the report will get bigger :)

![Overview Heatmap](./example_heatmap.png?raw=true "Overview heatmap showing quality criteria for each LC-MS file")
 
In the remainder of the report you can follow up on each individual metric and explore the reason for failure.

For example, the alignment via the Match-between-runs functionality is not optimal in this case, since MaxQuant cannot normalize the LC gradients at the start and end of the run:
![Alignment Performance](./example_MBRalignment.png?raw=true "Alignment of 5 raw files (the first file serves as reference here)")

Also, you could check how a new LC gradient influences the identifications over time:
![Identifications over Retention Time](./example_IDoverRT.png?raw=true "Identifications over Retention Time")

Find the full report here [Dataset2_5files_report_v0.61.0-100min.pdf].

See the package vignettes for documentation on how to create and customize a report.

  [Dataset2_5files_report_v0.61.0-100min.pdf]: Dataset2_5files_report_v0.61.0-100min.pdf

  
  [1]: http://www.ebi.ac.uk/pride/archive/projects/PXD003133


