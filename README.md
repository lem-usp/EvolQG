Morphometrics Package
======================

[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/lem-usp/Morphometrics?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Proto package.

Kinda usable


Instalation
===========

If Windows
----------

First install [Rtools](http://cran.r-project.org/bin/windows/Rtools/), then proceed to devtools instalation instructions.

Be careful to use the right version for your R instalation.


If MacOS X
-----------

If you are using MacOS X, you will probably have to build Rcpp from source. 
You will need to have [Xcode](https://developer.apple.com/xcode/) or some other C++ compiler installed.
Check if the instalation suceded typing `clang++ --version` in the terminal. It might be necessary to add a file called Makevars in a .R folder in your Home directory.

In the terminal type:

```
mkdir .R
cd .R
touch Makevars
```

edit this file (".R/Makevars") and add the following lines:

```
CC=clang
CXX=clang++
CFLAGS="-mtune=native -g -O3 -Wall -pedantic -Wconversion"
CXXFLAGS="-mtune=native -g -O3 -Wall -pedantic -Wconversion"
FLIBS=-lgfortran
```


Then, download the Rcpp source from [here](http://cran.r-project.org/web/packages/Rcpp/index.html) (Package source link), then install it in R using:

```R
> install.packages("./Rcpp_0.*.tar", repos = NULL, type="source")
```

Replacing the file name with the proper version number.

Then proceed to devtools instalation instructions.

Using devtools
--------------

Install [devtools](http://www.rstudio.com/projects/devtools/) and [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html), then:

```R
> library(devtools)
> install_github("lem-usp/Morphometrics")
```
