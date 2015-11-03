This folder contains the package vignettes in Markdown (.Rmd files).

When the package is build, R/knitr will automatically create 'proper' vignette files,
which end up as HTML in PTXQC/doc/*.html.

On GitHub, you can just open any of the vignettes and they will be displayed as HTML in your browser.


Technical thoughts on GitHub:

Since this package is also hosted on GitHub, vignettes can be directly viewed 'raw' by just browsing to this directory (PTXQC/vignettes).
GitHub has a nice feature, which renders the Markdown header as table (instead of displaying the YAML text or ugly LaTex commands).

Currently, GitHub will apply this nice formatting only if the vignette is encoded as ANSI. Using UTF-8 encoding showed the ugly header formatting.
Despite the fact that we promise UTF-8 to LaTex/Knitr in the header, ANSI seems to work nicely.

A markdown header is basically YAML, mixable with LaTex and might look like this:

---
title: "Basic R-Usage Guide for PTXQC"
author: "Chris Bielow <chris.bielow@mdc-berlin.de>"
date: '`r Sys.Date()`'
output:
  html_document: default
  pdf_document: null
vignette: >
  %\VignetteIndexEntry{Basic R-Usage Guide for PTXQC}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

See our vignettes for complete examples.