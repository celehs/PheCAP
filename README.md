
# PheCAP: High-Throughput Phenotyping with EHR using a Common Automated Pipeline

## Overview

The PheCAP package implements surrogate-assisted feature extraction
(SAFE) and common machine learning approaches to train and validate
phenotyping models. PheCAP begins with data from the EMR, including
structured data and information extracted from the narrative notes using
natural language processing (NLP). The standardized steps integrate
automated procedures, which reduce the level of manual input, and
machine learning approaches for algorithm training.

## Installation

The PheCAP package can be installed from CRAN or GitHub.

### Stable Version

Install stable version from CRAN:

``` r
install.packages("PheCAP")
```

### Development Version

Install development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("celehs/PheCAP")
```

## Get Started

Follow the [main
steps](https://celehs.github.io/PheCAP/articles/main.html), and try the
R codes from the [simulated
data](https://celehs.github.io/PheCAP/articles/example1.html) and [real
EHR data](https://celehs.github.io/PheCAP/articles/example2.html)
examples.

## References

-   Yichi Zhang`*`, Tianrun Cai`*`, Sheng Yu`*`, Kelly Cho, Chuan Hong,
    Jiehuan Sun, Jie Huang, Yuk-Lam Ho, Ashwin Ananthakrishnan, Zongqi
    Xia, Stanley Shaw, Vivian Gainer, Victor Castro, Nicholas Link,
    Jacqueline Honerlaw, Selena Huang, David Gagnon, Elizabeth Karlson,
    Robert Plenge, Peter Szolovits, Guergana Savova, Susanne Churchill,
    Christopher O’Donnell, Shawn Murphy, J Michael Gaziano, Isaac
    Kohane, Tianxi Cai`*`, and Katherine Liao`*`. [Methods for
    High-throughput Phenotyping with Electronic Medical Record Data
    Using a Common Semi-supervised Approach
    (PheCAP)](https://doi.org/10.1038/s41596-019-0227-6). *Nature
    Protocols* (2019). `*`contributed equally.

-   Yu, S., Chakrabortty, A., Liao, K. P., Cai, T., Ananthakrishnan, A.
    N., Gainer, V. S., … Cai, T. [Surrogate-assisted feature extraction
    for high-throughput
    phenotyping](https://doi.org/10.1093/jamia/ocw135). *Journal of the
    American Medical Informatics Association* (2017), e143-e149.

-   Liao, K. P., Cai, T., Savova, G. K., Murphy, S. N., Karlson, E. W.,
    Ananthakrishnan, A. N., … Kohane, I. [Development of phenotype
    algorithms using electronic medical records and incorporating
    natural language processing](https://doi.org/10.1136/bmj.h1885).
    *BMJ* (2015), 350(apr24 11), h1885–h1885.
