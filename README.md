PTXQC
---------------

[![License (3-Clause BSD)](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg?style=flat-square)](https://opensource.org/license/bsd-3-clause)
[![Project Stats](https://openhub.net/p/PTXQC/widgets/project_thin_badge.gif)](https://openhub.net/p/PTXQC)

**This package allows users of MaxQuant (from .txt files) and OpenMS (from mzTab files) to generate quality control reports in Html/PDF format.**

### Latest changes / ChangeLog

latest Release: v1.1.4 - Feb 2026<br>
latest Release on CRAN: same
  
See [NEWS][News_File] file for a version history.

### Platform support
  - Windows (recommended for convenience to make use of the [drag'n'drop][1] batch file provided)
  - Linux
  - MacOSX
  
### Citation

Please cite PTXQC when using it to check data in your publications:

**Proteomics Quality Control: Quality Control Software for MaxQuant Results**
Chris Bielow, Guido Mastrobuoni, and Stefan Kempa
*J. Proteome Res.*, 2016, 15 (3), pp 777-787.
DOI: [10.1021/acs.jproteome.5b00780][JPR_paper]
  
### Features
  - plethora of quality metrics
    - intensity distributions
    - digestion efficiency
    - contaminant visualizations
    - identification performance
    - Match-between-runs performance
  - easy usage ([Windows OS only] `drag'n'drop` your `txt output folder` onto a `batch file`)
    - 10 min [Installation](#installation)
  - MaxQuant locale check, i.e. dot as decimal separator (since PTXQC 1.0.10; see https://github.com/cbielow/PTXQC/pull/99 for details)
  - Html/PDF report will be generated within your MaxQuant-txt folder or next to the mzTab file
  - writes a mockup mzQC file (https://github.com/HUPO-PSI/mzQC/) for archiving or downstream processing (actual metrics still require exporting; PR's welcome)
  - optional configuration file *in YAML format* for generation of shorter/customized reports

### Target audience
  - MaxQuant users (no knowledge of R required)
  - OpenMS users (or any other software which can write an mzTab)
  - bioinformaticians (who want to contribute or customize)


### Documentation

A short **overview video** on PTX-QC] can be found [here](https://youtu.be/sb-mydbNRS4). 
If you want to jump to certain sections:

  - [Motivation for doing QC](https://youtu.be/sb-mydbNRS4?t=75)
  - [Installation demo](https://youtu.be/sb-mydbNRS4?t=522)
  - [Report Generation using Drag'n'Drop](https://youtu.be/sb-mydbNRS4?t=880)
  - [Inspecting and configuring a Report](https://youtu.be/sb-mydbNRS4?t=1043)
  
We use [pgkdown](https://pkgdown.r-lib.org/) to create a HTML documentation, which includes the Vignettes, our function documentation etc.
You can create the documentation locally: `pkgdown::build_site()` or visit the online-version on Github at [./docs/index.html](http://htmlpreview.github.io/?https://github.com/cbielow/PTXQC/blob/master/docs/index.html).

If you do not know where to start, look at the package vignettes first.

Our Vignettes give details on:

  - Full List of Quality Metrics with help text ('List of Metrics')
  - Input and Output
  - Report customization
  - (for MaxQuant/OpenMS users) Usage of Drag'n'drop
  - (for R users) Code examples in R

The 'List of Metrics' vignette contains a full description for each metric (the same as seen in the Help section of each Html report).

Within R, you can browse the vignettes using either of these commands (*after* the package is installed (see below)):
```
    help(package="PTXQC")
    browseVignettes(package = 'PTXQC')
```

Of course, you can also look at the Vignettes on CRAN:
[latest online vignette at CRAN](https://cran.r-project.org/package=PTXQC)

  
### Installation

**You can use PTX-QC without installing it, by using our webserver:**
Visit http://ptxqc.bsc.fu-berlin.de and simply upload your data.
This service is not suited for large-scale data analysis, but should be fine for the occasional analysis.

**If you want to generate QC reports without actually getting involved in R:**

We offer a Batch-file based Drag'n'drop mechanism to trigger PTXQC on any MaxQuant output folder.
This only works for Windows (not Linux or MacOS) at the moment -- but you have a Windows anyway to run MaxQuant, right?!
See [drag'n'drop][1] for details. It takes 10 minutes and you are done!

**If you just want to use PTXQC (and maybe even modify) it:**

*First*, install [pandoc](https://github.com/jgm/pandoc/releases) (see bottom of linked page). Pandoc is required in order to locally build the package vignettes (documentation),
but you can also read the [vignettes][Ref_Vign] online from the PTXQC GitHub page. More importantly, Pandoc enables PTXQC to write QC reports in HTML format (which come
with a help text for each plot and are interactive). PDF reports only contain plots!
The reports are printed as PDF by default and additionally as HTML if Pandoc is found.
**If you install Pandoc later while your R session is already open, you need to close and re-open R in order to make R aware of Pandoc!**
   
You can grab PTXQC from either [CRAN](https://CRAN.R-project.org/package=PTXQC) *or* [GitHub](https://github.com/cbielow/PTXQC#installation).
GitHub installation will give you the latest package; the CRAN version might be a little older, but is faster to install. Check the [NEWS][News_File] file for CRAN submissions and version.
> For the code blocks below: Run **each line** separately in your R console, i.e. do not copy and paste the whole block.
> If an error should occur, this allows to track it down more easily. See [FAQ - Installation][Ref_VignFAQ]
> how to resolve them.

    ## CRAN
    install.packages("PTXQC")
or

    ## GitHub
    if (!require(devtools, quietly = TRUE)) install.packages("devtools")
    library("devtools")             ## this might give a warning like 'WARNING: Rtools is required ...'. Ignore it.
    
    ## use build_vignettes = FALSE if you did not install pandoc or if you encounter errors when building vignettes (e.g. PRIDE ftp unavailable)!
    install_github("cbielow/PTXQC", build_vignettes = TRUE, dependencies = TRUE)

To get started, see the help and/or vignettes:

    help(package="PTXQC")
    browseVignettes(package = 'PTXQC')

Please feel free to report bugs (see below), or issue pull requests!    

### Report Examples

An overview chart at the beginning of the report will give you a first impression.
<img src="./inst/examples/example_heatmap.png?raw=true" width="500" /><br>
Detailed plots can be found in the remainder of each report.

For example input data and full reports, see the ['inst/examples' subfolder][2].

### Bug reporting / Feature requests

If you encounter a bug, e.g. error message, wrong figures, missing axis annotation or anything which looks
suspicious, please use the [GitHub issue tracker][issuetracker] and file a report.

You should include
  - **stage** you encounter the bug, e.g. during installation, report creation, or after report creation (i.e. a bug in the report itself).
  - **PDF/Html report** itself (if one was generated).
  - **version of PTXQC**, e.g. see the report_XXX.pdf/html (where XXX will be the version) or see the DESCRIPTION file of the PTXQC package or call `help(package="PTXQC")` within R
  - **error message** (very important!). Either copy it or provide a screen shot.

Please be as precise as possible when providing the bug report - just imagine what kind of information you would like to have in order
to track down the issue.
In certain situations, the whole txt-folder or a single MaxQuant/mzTab file might be helpful to solve the problem.

### Contributing - Get Involved!

We welcome input from our user base!
PTX-QC has a very permissive **BSD-3 clause License** (see [DESCRIPTION](DESCRIPTION) file), so feel free to fork, patch and contribute!
There are many ways to get involved, _you do not need to be a developer_!
  - suggest a new metric (and why you think it's useful) by opening [a new ticket][issuetracker] here on GitHub.
  - suggest changes to existing metrics (improvements or bugfixes), see above.
  - suggest improvements to our documentation (e.g. [additional vignettes][Ref_Vign])
  - write code (in R) and submit a [Pull Request (PR)][PullRequest].

### PTX-QC user forum

Come and interact with developers and users at https://discord.gg/MB6PvpctUY


### Misc

Use [PTXQC v0.69.3][JPR_PTXQC] if you want the version which was used in the paper, i.e.
use `install_github(..., ref="v0.69.3")` when following the [Installation](#installation) procedure.

The input data used in the original publication is available in the ['inst/examples' subfolder][2].

We recommend to use the most recent PTXQC for the best user experience.
  
  [1]: https://github.com/cbielow/PTXQC/tree/master/inst/dragNdrop
  [2]: https://github.com/cbielow/PTXQC/tree/master/inst/examples
  [issuetracker]: https://github.com/cbielow/PTXQC/issues
  [PullRequest]: https://github.com/cbielow/PTXQC/pulls
  [JPR_PTXQC]: https://github.com/cbielow/PTXQC/releases/tag/v0.69.3
  [Ref_VignFAQ]: https://github.com/cbielow/PTXQC/blob/master/vignettes/PTXQC-FAQ.Rmd
  [Ref_Vign]: https://github.com/cbielow/PTXQC/tree/master/vignettes
  [News_File]: https://github.com/cbielow/PTXQC/blob/master/NEWS
  [JPR_paper]: https://doi.org/10.1021/acs.jproteome.5b00780
