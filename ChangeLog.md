# Change Log
All notable changes to `mptp` will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.2.5] - 2023-09-11
### Added
 - Likelihood ratio test for the multi method
 - Added implementation for the incomplete gamma function
### Removed
 - Dependency for GNU scientific library

## [0.2.4] - 2018-05-14
### Fixed
 - If we do not manage to generate a random starting delimitation with the
   wanted number of species (randomly chosen), we use the currently generated
   delimitation instead.

## [0.2.3] - 2017-07-25
### Fixed
 - Replaced hsearch which was causing problems on APPLE with custom hashtable
 - Corrected file name in error messages when failing to open files

## [0.2.2] - 2017-01-31
### Fixed
 - Regular expressions now allow scientific notation when parsing branch lengths
 - Improved accuracy of ASV score (takes into account tip species)
 - Memory leaks when parsing incorrectly formatted trees

## [0.2.1] - 2016-10-18
### Fixed
 - Updated ASV to consider only coalescent roots of ML delimitation
 - Assertion stopping mptp when using random starting delimitations for MCMC

## [0.2.0] - 2016-09-27
### Fixed
 - Floating point exception error when constructing random trees caused from
   division by zero
 - Allocation with malloc caused uninitialized variables when converting unrooted
   tree to rooted for the MCMC method
 - Sample size for the the AIC with a correction for finite sample sizes
