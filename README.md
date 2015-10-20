PTXQC
---------------

[![Project Stats](https://www.openhub.net/p/PTXQC/widgets/project_thin_badge.gif)](https://www.openhub.net/p/PTXQC)

**This package allows users of MaxQuant to generate quality control reports in PDF format.**

### Features
  - plethora of quality metrics
    - intensity distributions
    - digestion efficiency
    - contaminant visualizations
    - identification performance
    - Match-between-runs performance
  - easy usage ([Windows OS only] `drag'n'drop` your `txt output folder` onto a `batch file`)
    - 10 min installation, see [drag'n'drop] [1]
  - PDF report will be generated within your txt folder
  - optional configuration file *in YAML format* for generation of shorter/customized reports

### Target audience
  - MaxQuant users (no knowledge of R required)
  - bioinformaticians (who want to contribute or customize)

### Installation

**If you want to generate QC reports without actually getting involved in R:**

We offer a Batch-file based Drag'n'drop mechanism to trigger PTXQC on any MaxQuant output folder.
This only works for Windows (not Linux or MacOS) at the moment -- but you have a Windows anyway to run MaxQuant, right?!
See [drag'n'drop] [1] for details. It takes 10 minutes and you are done!

**If you just want the package to use and maybe even modify it:**

Direct installation from GitHub requires the 'devtools' package. The following lines will add PTXQC on a fresh R installation:

    if (!require(devtools, quietly = TRUE)) install.packages("devtools")
    library("devtools")             ## this might give a warning like 'WARNING: Rtools is required ...'. Ignore it.
    source("http://bioconductor.org/biocLite.R")
    biocLite("Biobase")
    install_github("cbielow/PTXQC", build_vignettes = TRUE) 
    
To get started on how to use the package, see the help and/or vignettes:

    help(package="PTXQC")
    browseVignettes(package = 'PTXQC')

Please feel free to report bugs (see below), or issue pull requests!    
    
### Platform support

  - Windows (recommended for convenience to make use of the [drag'n'drop] [1] batch file provided)
  - Linux
  - MacOSX

### Bug reporting / Feature requests

If you encounter a bug, e.g. error message, wrong figures, missing axis annotation or anything which looks
suspicious, please use the GitHub issue tracker (top right) and file a report.

You should mention the **version of PTXQC** which contains the bug (see the report_XXX.pdf, where XXX will be the version),
and the **PDF report** itself (if one was generated).
Also the **error message** you see is probably very important. Either copy it or provide a screen shot.
Please be a precise as possible when providing the bug report -- just imagine what kind of information you would like to have in order
to track down the issue.
In certain situations, the whole txt-folder or a single MaxQuant file might be helpful to solve the problem.

If you want to see a new metric, or have ideas how to improve the existing ones, just open an issue ticket and leave a description.
  
### Report Examples

An overview chart at the beginning of the report will give you a first impression.
<img src="./inst/examples/example_heatmap.png?raw=true" width="500" /><br>
Detailed plots can be found in the remainder of each report.

For more, see the ['inst/examples' subfolder] [2].

  
  [1]: https://github.com/cbielow/PTXQC/tree/master/inst/dragNdrop
  [2]: https://github.com/cbielow/PTXQC/tree/master/inst/examples
