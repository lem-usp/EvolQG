Morphometrics Package
======================

Proto package.

Kinda usable


Instalation
===========

If Windows
----------

Install [Rtools](http://cran.r-project.org/bin/windows/Rtools/). 

Be careful to use the right version for your R instalation.


If MacOS X
-----------

If you are using MacOS X, you will probably have to build Rcpp from source. 
You will need to have [Xcode](https://developer.apple.com/xcode/) or some other C++ compiler installed.
Check if the instalation suceded typing `clang++ --version` in the terminal.


Then, download the Rcpp source from [here](http://cran.r-project.org/web/packages/Rcpp/index.html) (Package source link), then install it in R using:

```R
> install.packages("Rcpp_0.*.tar.gz", repos = NULL, type="source")
```

Replacing the file name with the proper version number.


Using devtools
--------------

Install [devtools](http://www.rstudio.com/projects/devtools/) and [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html), then:

```R
> library(devtools)
> install_github("Morphometrics", "lem-usp")
```

Deixando o Morphometrics como padrão no RStudio
-----------------------------------------------

1. localizar o arquivo Rprofile.site. No windows ele fica no C:\Program Files\R\R-n.n.n\etc 
2. abrir o RStudio como administrador
3. abrir o Rprofile.site no RStudio
4. após a última linha com texto copiar options(defaultPackages=c(getOption("defaultPackages"),"Morphometrics"))
5. salvar e fechar o Rprofile.site
6. reiniciar o RStudio
7. pronto!
