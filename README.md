# mvregrp

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**Multivariate Grouped Random Effect Model Estimation**

The `mvregrp` package implements efficient Gibbs sampling for multivariate grouped random effect models, designed for analyzing household panel data with time-varying household membership and correlation between household random effects within 'superhousehold' clusters.

## Key Features

- **Grouped Random Effects**: Household effects that adapt to changing composition over time
- **Correlated Random Effects**: Models correlations between household effects within superhousehold clusters  
- **Covariate-Dependent Random Effect Correlations**: Correlation structure depends on household relationship covariates
- **Robust MCMC**: Novel algorithms ensuring positive definite correlation matrices
- **Multiple Variants**: Models with/without area effects and different household structures

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("slzhang-fd/mvregrp")
```

## Test Data

The package includes a test dataset (`test.rda`) that demonstrates the model functionality. The data were generated following the design of the simulation study described in Section 6.1 with heterogeneous between-household correlation structure. The dataset contains 1000 superhouseholds and includes the following components:

- **ysim**: Binary outcome variable
- **x_covs**: Mean model design matrix
- **corr_covs**: Correlation model design matrix 
- **i_ind**: Individual indices
- **sh_ind**: Super-household indices
- **hit_ind**: Household indices

### Data Structure Examples

The `x_covs` variables correspond to x0-x4 in Section 6.1 of the paper:
```
head(x_covs):
     cons age50d10 highed tenurehh2 tenurehh3
[1,]    1      1.0      0         0         0
[2,]    1      1.1      0         0         0
[3,]    1      1.2      0         0         0
[4,]    1      1.3      0         0         0
[5,]    1      1.2      0         0         0
[6,]    1      1.4      0         0         0
```

The `corr_covs` contain correlation model covariates. The first three variables are the superhousehold ID and the household IDs for each pair of observations in a superhousehold (where there is more than one household, and the order of corr_covs records does not matter). The other variables correspond to z0-z4 in Section 6.1 of the paper:
```
head(corr_covs):
  sh_ind hit_ind1 hit_ind2 cons part1parch0 part0parch1 part0parch0      pbothc
1      2     1377     1378    1           1           0           0  0.00000000
2      4     1579     1580    1           1           0           0  0.30000001
3      4     1579     1581    1           0           0           0  0.10000002
4      4     1579     1582    1           0           0           0 -0.35714285
5      4     1579     1583    1           0           0           0 -0.09999999
6      4     1579     1584    1           0           1           0 -0.25000000
```

## Usage

Load example data and fit a grouped random effect model:

```r
library(mvregrp)
library(matrixStats)  # For colSds, colQuantiles
library(coda)         # For effectiveSize

# Load the test dataset
data_path <- system.file("extdata", "test.rda", package = "mvregrp")
load(data_path)

# Set MCMC parameters
nsample <- 10000
burn_in <- 5000

# Run grouped random effect model
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

Note that all specifications include individual random effects.

| Function | Description |
|----------|-------------|
| `svregrp_Gibbs()` | Standard grouped household effect model |
| `svregrp_Gibbs_area()` | Includes additional area random effects |
| `svregrp_Gibbs_area_nohe()` | Area effects only (no household effects) |
| `svregrp_Gibbs_mchains()` | Multiple chains for convergence diagnostics |

## Model Framework

The package implements random effects models for longitudinal data with the following complex non-hierarchical structure:

- **Time-varying household composition**: Household membership can change over time  
- **Correlated random effects**: Household random effects correlated within superhousehold clusters
- **Flexible between-household correlations**: Correlation structure depends on household relationship covariates

Details of the data structure are given after eq. (7) of the paper.

Applications include health outcomes, social behaviors, and economic variables in studies like the UK Household Longitudinal Study (UKHLS).

## Citation

Please cite this package as:

> Steele, F., Zhang, S., & Clarke, P. (2025). Analysis of household effects on longitudinal health outcomes using a joint mean-correlation multilevel model with grouped random effects. *Manuscript under review*.

## License

GPL-3
