[![cran version](http://www.r-pkg.org/badges/version/baorista)](https://CRAN.R-project.org/package=baorista) 
[![development version](https://img.shields.io/badge/devel%20version-0.1.2-lightblue.svg)](https://github.com/ercrema/baorista)
[![R-hub](https://github.com/ercrema/baorista/actions/workflows/rhub.yaml/badge.svg)](https://github.com/ercrema/baorista/actions/workflows/rhub.yaml)

# baorista  <img src="/logo/logo.png" align="right" />
_baorista_ is an R package that provides a Bayesian inferential tool for analysing time-frequencies of archaeological events associated with time spans typically obtained from relative chronological sequences (e.g. periods and phases) and often analyses using aoristic sums. At its core _baorista_ is a frontend for fitting Bayesian models via the [NIMBLE probabilistic programming language](https://r-nimble.org/). The package is currently in beta development and can be installed using `devtools`:

```
library(devtools)
install_github('ercrema/baorista')
```

Please note that baorista is based on Nimble, which requires a working C++ compiler. For more information please read the [dedicated section of the nimble manual](https://r-nimble.org/html_manual/cha-installing-nimble.html#sec:compiler).

For a quick introduction to _baorista_ check the package [vignette](https://htmlpreview.github.io/?https://github.com/ercrema/baorista/blob/main/vignettes/using_baorista.html) as well as the associated [paper](https://doi.org/10.1111/arcm.12984) and its [github repo](https://github.com/ercrema/beyond_aoristic)

## Funding
The development of this package was funded by a Philip Leverhulme Prize (PLP-2019-304) awarded to E.Crema.

