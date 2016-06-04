EvolQG Package
======================

[![Join the chat at https://gitter.im/diogro/evolqg](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/diogro/evolqg?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Travis build status](https://travis-ci.org/lem-usp/EvolQG.svg?branch=master)](https://travis-ci.org/lem-usp/EvolQG) 
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/lem-usp/EvolQG?branch=master&svg=true)](https://ci.appveyor.com/project/lem-usp/EvolQG)
[![Coverage Status](https://coveralls.io/repos/github/lem-usp/EvolQG/badge.svg?branch=master)](https://coveralls.io/github/lem-usp/EvolQG?branch=master)

Package for evolutionary quantitative genetics, described in Melo D, Garcia G, Hubbe A et al. EvolQG - An R package for evolutionary quantitative genetics [version 1; referees: 1 approved, 1 approved with reservations] F1000Research 2015, 4:925 (doi: 10.12688/f1000research.7082.1) [article link](http://f1000research.com/articles/4-925/v1)


Installation
============

From CRAN
---------

```R
> install.packages("evolqg")
```

From github
-----------

In Windows you first need to install [Rtools](http://cran.r-project.org/bin/windows/Rtools/), then proceed to devtools installation instructions. Be careful to use the right version for your R installation.


Using devtools
--------------

Install [devtools](http://www.rstudio.com/projects/devtools/), then:

```R
> library(devtools)
> install_github("lem-usp/evolqg")
```

Troubleshooting
---------------

All dependencies should be available from CRAN. Please open an [issue](https://github.com/lem-usp/EvolQG/issues) if install fails for some reason.
