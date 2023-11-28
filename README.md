### GeneBridge: Rooting Orthologous Genes in Large-Scale Evolutionary Analysis.

*GeneBridge* is an *R* package designed to assess the evolutionary history of genes across diverse species. It implements the Bridge algorithm to infer the evolutionary root of genes in a given species tree. By rooting orthologous genes in large-scale evolutionary snalysis, *GeneBridge* can provide a framework for exploring genes within biological systems.

### Installation in R (>=4.3)

##### Install dependencies to build the package's vignettes

```r
install.packages("knitr")
install.packages("rmarkdown")
install.packages("BiocManager")
BiocManager::install("BiocStyle")
```

##### Install the GeneBridge package

```r
install.packages("remotes")
remotes::install_github("sysbiolab/GeneBridge", build_vignettes=TRUE)
```

### Examples

Follow the *GeneBridge* vignette and try to make some *plots*!

```r
library(GeneBridge)
vignette("GeneBridge")
```

### Licenses

The *GeneBridge* package is distributed under [Artistic-2.0](https://www.r-project.org/Licenses/Artistic-2.0)
