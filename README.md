
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multimetab

<!-- badges: start -->

<!-- badges: end -->

`multimetab` is an R package for performing Bayesian variable selection
with mass spectrometry-based metabolite measurements as the response
variable. It provides companion code for our manuscript, “A
nutritionally informed model for Bayesian variable selection of
metabolite endpoints,” which presented a novel framework metabolite
analysis based on a skew-normal censored mixture model (SNCM).
Metabolite measurements require specialized models because they are
highly right-skewed and exhibit missing values, called point mass values
(PMVs), that arise from multiple sources. “Biological” PMVs occur if the
mass spectrometer fails to detect the metabolite because it is
completely absent from the sample in question. “Technical” PMVs occur
when the metabolite is present in the sample but with an intensity below
the metabolite’s “detection limit”, resulting in a spurious missing
value that does not reflect total absence. The proposed internally
handles PMVs of both types, while directly modeling the skewness of the
positive-valued metabolites via the skew-normal distribution. Variable
selection is performed by placing spike-and-slab priors on the
regression coefficients and, optionally, using a Markov Random Field
(MRF) prior that exploits external information on the substantive
relationships among variables.

## Installation

You can install the development version of multimetab from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dclarkboucher/multimetab")
```
