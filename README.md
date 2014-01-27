Morphometrics Package
======================

Proto package.

Kinda usable


Instalation
===========

If Windows
----------

Install [Rtools](http://cran.r-project.org/bin/windows/Rtools/)


Using devtools
--------------

Install [devtools](http://www.rstudio.com/projects/devtools/), then:

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
