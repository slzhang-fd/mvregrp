# mvregrp 1.0.0

## Initial Release

* Initial release of mvregrp package
* Implements Gibbs sampling for multivariate grouped random effect probit models
* Supports time-varying household membership with correlated random effects
* Includes three main modeling functions:
  - `svregrp_Gibbs()`: Standard grouped household effect model
  - `svregrp_Gibbs_area()`: Model with area random effects
  - `svregrp_Gibbs_area_nohe()`: Simplified model without household effects
* Novel MCMC estimation ensuring positive definite correlation matrices
* Covariate-dependent correlations between household effects
* Multiple chain support for convergence diagnostics
* Utility functions for model comparison and visualization