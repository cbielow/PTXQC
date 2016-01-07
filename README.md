PTXQC
---------------

[![Build Status](https://travis-ci.org/cbielow/PTXQC.svg?branch=master)](https://travis-ci.org/cbielow/PTXQC) 
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

### Documentation
  
Besides this documentation on GitHub, the package vignettes
of PTXQC will give you valuable information. *After* the package is installed (see below),
you can browse the vignettes using either of these commands within R:

    help(package="PTXQC")
    browseVignettes(package = 'PTXQC')
  
If you do not want to wait that long, have a look at the ['vignettes' subfolder] [3].
The top part contains a small table with technical gibberish, but the rest is identical to the
vignettes you would see in R.

You will find documentation on
  - Input and Output
  - Report customization
  - (for MaxQuant users) Usage of Drag'n'drop
  - (for R users) code examples in R
  - ...
  
### Installation

**If you want to generate QC reports without actually getting involved in R:**

We offer a Batch-file based Drag'n'drop mechanism to trigger PTXQC on any MaxQuant output folder.
This only works for Windows (not Linux or MacOS) at the moment -- but you have a Windows anyway to run MaxQuant, right?!
See [drag'n'drop] [1] for details. It takes 10 minutes and you are done!

**If you just want the package to use (and maybe even modify) it:**

The following lines will install the PTXQC package.
Direct installation from GitHub requires the 'devtools' package.
Run **each line** separately in your R console, i.e. do not copy and paste the whole block.
If an error should occur, this allows to track it down more easily. See [FAQ - Installation] [Ref_VignFAQ]
how to resolve them.
   
    if (!require(devtools, quietly = TRUE)) install.packages("devtools")
    library("devtools")             ## this might give a warning like 'WARNING: Rtools is required ...'. Ignore it.
    source("http://bioconductor.org/biocLite.R")
    biocLite("Biobase")
    install_github("cbielow/PTXQC", build_vignettes = TRUE) 
    
To get started, see the help and/or vignettes:

    help(package="PTXQC")
    browseVignettes(package = 'PTXQC')

Please feel free to report bugs (see below), or issue pull requests!    
    
### Platform support

  - Windows (recommended for convenience to make use of the [drag'n'drop] [1] batch file provided)
  - Linux
  - MacOSX

### Bug reporting / Feature requests

If you encounter a bug, e.g. error message, wrong figures, missing axis annotation or anything which looks
suspicious, please use the [GitHub issue tracker] [issuetracker] and file a report.

You should include
  - **stage** you encounter the bug, e.g. during installation, report creation, or after report creation (i.e. a bug in the report itself).
  - **PDF report** itself (if one was generated).
  - **version of PTXQC**, e.g. see the report_XXX.pdf (where XXX will be the version) or see the DESCRIPTION file of the PTXQC package or call `help(package="PTXQC")` within R
  - **error message** (very important!). Either copy it or provide a screen shot.

Please be as precise as possible when providing the bug report - just imagine what kind of information you would like to have in order
to track down the issue.
In certain situations, the whole txt-folder or a single MaxQuant file might be helpful to solve the problem.

If you want to see a new metric, or have ideas how to improve the existing ones, just open an issue ticket and leave a description.
  
### Report Examples

An overview chart at the beginning of the report will give you a first impression.
<img src="./inst/examples/example_heatmap.png?raw=true" width="500" /><br>
Detailed plots can be found in the remainder of each report.

For example input data and full reports, see the ['inst/examples' subfolder] [2].

  
### Citation

PTXQC is published at JPR:
**Proteomics Quality Control: Quality Control Software for MaxQuant Results**
Chris Bielow, Guido Mastrobuoni, and Stefan Kempa
*J. Proteome Res.*
Publication Date (Web): December 14, 2015
DOI: [https://doi.org/10.1021/acs.jproteome.5b00780](10.1021/acs.jproteome.5b00780)

Use [PTXQC v0.69.3] [JPR_PTXQC] if you want the version which was used in the paper, i.e.
use `install_github("cbielow/PTXQC", ref="v0.69.3", build_vignettes = TRUE)` when following the [Installation](#Installation) procedure.

The input data is available in the ['inst/examples' subfolder] [2].

We recommend to use the most recent PTXQC for the best user experience.

  
  [1]: https://github.com/cbielow/PTXQC/tree/master/inst/dragNdrop
  [2]: https://github.com/cbielow/PTXQC/tree/master/inst/examples
  [3]: https://github.com/cbielow/PTXQC/tree/master/vignettes
  [issuetracker]: https://github.com/cbielow/PTXQC/issues
  [JPR_PTXQC]: https://github.com/cbielow/PTXQC/releases/tag/v0.69.3
  [Ref_VignFAQ]: https://github.com/cbielow/PTXQC/blob/master/vignettes/PTXQC-FAQ.Rmd