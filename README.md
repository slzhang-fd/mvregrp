# mvregrp: Multivariate Grouped Random Effect Probit Model Estimation

<!-- badges: start -->
<!-- badges: end -->

The `mvregrp` package implements efficient Gibbs sampling procedures for multivariate grouped random effect probit model estimation. This package is designed for analyzing longitudinal data with time-varying household membership, allowing for correlated random effects within 'superhousehold' clusters.

## Key Features

- **Grouped Random Effects**: Models household effects that can change over time as household composition changes
- **Correlated Effects**: Allows correlation between household random effects within superhousehold clusters
- **Covariate-Dependent Correlations**: Models correlations as functions of covariates describing household relationships
- **Robust MCMC**: Novel MCMC estimation procedures that ensure positive definiteness of correlation matrices
- **Multiple Model Variants**: Supports models with and without area effects, and with different household structures


## Model Variants

The package provides several modeling functions:

- `svregrp_Gibbs()`: Standard grouped household effect model
- `svregrp_Gibbs_area()`: Model with additional area random effects  
- `svregrp_Gibbs_area_nohe()`: Simplified model without household effects
- `svregrp_Gibbs_mchains()`: Multiple chain implementation for convergence diagnostics

## Model Description

The package implements a hierarchical model for longitudinal data where:

- Individuals are nested within households
- Households can change composition over time
- Related households are grouped into "superhouseholds" 
- Random effects are allowed to be correlated within superhousehold clusters
- Correlations can depend on covariates describing household relationships

This is particularly useful for analyzing health outcomes, social behaviors, or economic variables in household panel studies like the UK Household Longitudinal Study (UKHLS).


## Installation & Usage

You can install the development version of mvregrp from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("slzhang-fd/mvregrp")
```
### Basic Example

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

## Citation

If you use this package in your research, please cite:

> Steele, F., Zhang, S., & Clarke, P. (2025). Analysis of household effects on longitudinal health outcomes using a joint mean-correlation multilevel model with grouped random effects. *Manuscript under review*.

## License

GPL-3