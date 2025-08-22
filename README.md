# mvregrp

<!-- badges: start -->
[![R-CMD-check](https://github.com/slzhang-fd/mvregrp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/slzhang-fd/mvregrp/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**Multivariate Grouped Random Effect Model Estimation**

The `mvregrp` package implements efficient Gibbs sampling for multivariate grouped random effect models, designed for analyzing longitudinal household panel data with time-varying membership and correlated effects within clustered 'superhouseholds'.

## Key Features

- **Grouped Random Effects**: Household effects that adapt to changing composition over time
- **Correlated Effects**: Models correlations between household effects within superhousehold clusters  
- **Covariate-Dependent Correlations**: Correlation structure depends on household relationship covariates
- **Robust MCMC**: Novel algorithms ensuring positive definite correlation matrices
- **Multiple Variants**: Models with/without area effects and different household structures

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("slzhang-fd/mvregrp")
```

## Usage

Load example data and fit a grouped random effect model:

```r
library(mvregrp)
library(matrixStats)  # For colSds, colQuantiles
library(coda)         # For effectiveSize

# Load pre-prepared test data
# This includes: ysim (outcome), x_covs (mean model covariates), 
# corr_covs (correlation model covariates), and index vectors
data_path <- system.file("extdata", "test.rda", package = "mvregrp")
load(data_path)

# Set MCMC parameters
nsample <- 10000
burn_in <- 5000

# Run grouped random effect probit model
result <- svregrp_Gibbs(
  Y_star = as.matrix(ysim),           # Outcome variable
  x_covs = as.matrix(x_covs),         # Mean model design matrix
  z_covs = as.matrix(corr_covs),      # Correlation model design matrix
  i_ind = i_ind,                      # Individual indices
  sh_ind = sh_ind,                    # Super-household indices  
  hit_ind = hit_ind,                  # Household indices
  max_steps = nsample,                # Total MCMC iterations
  cor_step_size = c(rep(0.01, ncol(corr_covs) - 3), 0.001),  # MH jump step sizes for correlation parameters, needs pre-tunning
  corr_vs_diag = FALSE,               # Use full correlation structure
  verbose = FALSE                     # Suppress progress output
)

# Extract and process MCMC chains
chains <- as.data.frame(result$params_mcmc_obj)

# Reorder columns: mean coefficients, variance parameters, correlation coefficients
chains_ordered <- as.matrix(cbind(
  chains[1:5],        # Mean model coefficients
  chains[11:13],      # Variance parameters (sigma2_e, sigma2_u, sigma2_v)
  chains[6:10]        # Correlation coefficients
))

# Compute posterior summaries (excluding burn-in)
post_samples <- chains_ordered[(burn_in + 1):nsample, ]
posterior_means <- colMeans(post_samples)
posterior_sd <- colSds(post_samples)
credible_intervals <- colQuantiles(post_samples, probs = c(0.025, 0.975))
effective_sample_size <- effectiveSize(post_samples)

# Display results
print("Posterior means:")
print(round(posterior_means, 3))

print("95% Credible intervals:")
print(round(credible_intervals, 3))

print("Effective sample sizes:")
print(round(effective_sample_size))
```

## Model Variants

| Function | Description |
|----------|-------------|
| `svregrp_Gibbs()` | Standard grouped household effect model |
| `svregrp_Gibbs_area()` | Includes additional area random effects |
| `svregrp_Gibbs_area_nohe()` | Area effects only (no household effects) |
| `svregrp_Gibbs_mchains()` | Multiple chains for convergence diagnostics |

## Model Framework

The package implements hierarchical models for longitudinal household data with:

- **Nested structure**: Individuals within households within superhouseholds
- **Time-varying composition**: Household membership can change over time  
- **Correlated effects**: Household random effects correlated within superhousehold clusters
- **Flexible correlations**: Correlation structure depends on household relationship covariates

Applications include health outcomes, social behaviors, and economic variables in studies like the UK Household Longitudinal Study (UKHLS).

## Citation

Please cite this package as:

> Steele, F., Zhang, S., & Clarke, P. (2025). Analysis of household effects on longitudinal health outcomes using a joint mean-correlation multilevel model with grouped random effects. *Manuscript under review*.

## License

GPL-3
