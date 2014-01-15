Morphometrics Package
======================

Proto package.

Kinda usable


Instalation
===========

Using devtools
--------------

Install [devtools](http://www.rstudio.com/projects/devtools/), then:

```R
> library(devtools)
> install_github("Morphometrics", "lem-usp")
```


Using package tar.gz
--------------------

Download [this file](https://dl.dropboxusercontent.com/u/891794/Morphometrics_0.1.tar.gz) and run:

```R
> install.packages(c("plyr", "reshape2", "ggplot2", "vegan", "mvtnorm"))
> install.packages("./Morphometrics_0.1.tar.gz")
> library(Morphometrics)
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
