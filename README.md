PTXQC
---------------

**This package allows users of MaxQuant to generate quality control reports in PDF format.**

### Features
  - plethora of quality metrics
    - intensity distributions
    - digestion efficiency
    - contaminant visualizations
  - easy usage (drag'n'drop your txt output folder onto a batch file)
    - 10 min installation, see ['inst' subfolder] [1] 
  - PDF report will be generated within your txt folder
  - optional configuration file for generation of shorter/customized reports

### Target audience
  - MaxQuant users (no knowledge of R required)
  - bioinformaticians (who want to contribute or customize)

### Installation
Direct installation from GitHub requires the 'devtools' package. This should install PTXQC on a fresh R system

    install.packages("devtools")
    library("devtools")
    install_github("cbielow/PTXQC")


### Platform support
  - Windows (recommended for convenience to make use of the drag'n'drop batch file provided)
  - Linux
  - MacOSX

  
  [1]: https://github.com/cbielow/PTXQC/inst/dragNdrop/
