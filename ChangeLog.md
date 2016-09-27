# Change Log
All notable changes to `mptp` will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.2.0] - 2016-09-27
### Fixed
 - Floating point exception error when constructing random trees caused from
   division by zero
 - Allocation with malloc caused uninitialized variables when converting unrooted
   tree to rooted for the MCMC method
 - Sample size for the the AIC with a correction for finite sample sizes
