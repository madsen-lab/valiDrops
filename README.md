# valiDrops

valiDrops is an R package for detecting high-quality barcodes in droplet-based single-cell or single-nucleus RNA-seq datasets. In addition, valiDrops is capable of detecting and flagging apoptotic cells.

# Installation

valiDrops can be installed directly from GitHub using either {devtools} or {remotes}. In addition to valiDrops, it is high advisable to install {presto}, which yields a 80 - 100X speed-up when using valiDrops. Both packages can be installed using the following commands:

```{r}
install.packages("remotes")
remotes::install_github("madsen-lab/valiDrops")
remotes::install_github("immunogenomics/presto")
```

# Usage
