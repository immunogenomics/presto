# Overview

Presto performs a fast Wilcoxon rank sum test and auROC analysis. Latest benchmark ran 1 million observations, 1K features, and 10 groups in 16 seconds (sparse input) and 85 seconds (dense input). 


# Installation

We are working on getting Presto into CRAN. For now, install Presto from github directly:

```{r}
library(devtools)
install_github('immunogenomics/presto')
```

# Usage

Run presto on a matrix, Seurat, or SingleCellExperiment input object. 

```
wilcoxauc(X, y)
wilcoxauc(seurat_object, 'group_name')
wilcoxauc(sce_object, 'group_name')
```

For examples, see `?wilcoxauc` and the [vignette](http://htmlpreview.github.io/?https://github.com/immunogenomics/presto/blob/master/docs/getting-started.html)
