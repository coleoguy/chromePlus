
# chromePlus

<!-- badges: start -->
[![R-CMD-check](https://github.com/coleoguy/chromePlus/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/coleoguy/chromePlus/actions/workflows/R-CMD-check.yaml)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
<!-- badges: end -->

An R package for modeling chromosome number evolution with binary
trait-dependent rates. chromePlus extends the chromEvol framework by
incorporating a binary trait that can affect rates of chromosome gain, loss,
polyploidy, and demiploidy, as well as state-dependent speciation and
extinction. Built on the [diversitree](https://cran.r-project.org/package=diversitree)
framework.

## Installation

```r
# Install from GitHub
install.packages("remotes")
remotes::install_github("coleoguy/chromePlus")
```

## Quick example

```r
library(chromePlus)
library(diversitree)

# Simulate chromosome evolution on a random tree
set.seed(123)
tree <- rcoal(30)
result <- simChrom(
  tree = tree,
  pars = c(0.1, 0.08, 0.0, 0.01, 7),
  limits = c(3, 15),
  model = "2010"
)
head(result)
```

## Key features

- **Multiple evolution models**: Simple chromEvol (2010), ChromPlus with binary
  traits, PloidEvol with ploidy as a hidden state, and sex chromosome-autosome
  fusion (SAF) model
- **Flexible constraints**: Drop polyploidy/demiploidy, force symmetric rates,
  constrain hyperstates, and more
- **State-dependent diversification**: MuSSE integration allows testing whether
  chromosome evolution correlates with speciation/extinction rates
- **Simulation tools**: `simChrom()` for fixed-tree simulations and
  `makeSSEchrom()` for full SSE simulations
- **MCMC visualization**: `plotChromeplus()` for density plots with HPD intervals
- **Uncertainty support**: Handle uncertain binary trait assignments via
  probability matrices

## Citation

To cite chromePlus in publications, use:

```r
citation("chromePlus")
```

Blackmon, H., Justison, J., Mayrose, I. and Goldberg, E.E. (2019). Meiotic
drive shapes rates of karyotype evolution in mammals. *Evolution*, 73(3),
511--523.

Blackmon, H., Chin, M. and Jonika, M.M. (2023). chromePlus: Analysis of
Chromosome Number Evolution and Binary Traits.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8102439.svg)](https://doi.org/10.5281/zenodo.8102439)

## Issues

If you have questions or encounter problems, please open an issue at
<https://github.com/coleoguy/chromePlus/issues>.
