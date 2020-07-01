---
layout: page
title: Installation
---

For installation of the [CellNOptR package](https://bioconductor.org/packages/release/bioc/html/CellNOptR.html), open `R` and type:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CellNOptR", version = "3.8")
```


**Other CNO formalisms are extension of the CellNOptR package and can only be installed after CellNOptR is present on the library of `R` packages:**


To install [CNORdt package](https://bioconductor.org/packages/release/bioc/html/CNORdt.html), start `R` and enter:

```R
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CNORdt", version = "3.8")
```



To install [CNORfeeder package](https://bioconductor.org/packages/release/bioc/html/CNORfeeder.html), start `R` and enter:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CNORfeeder", version = "3.8")
```



To install [CNORfuzzy package](https://bioconductor.org/packages/release/bioc/html/CNORfuzzy.html), start `R` and enter:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CNORfuzzy", version = "3.8")
```


To install [CNORode package](https://bioconductor.org/packages/release/bioc/html/CNORode.html), start `R` and enter:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CNORode", version = "3.8")
```
