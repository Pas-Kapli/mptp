# Species Delimitation

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)
[![Build Status](https://travis-ci.com/Pas-Kapli/mptp.svg?token=Rni6Dxb79oRVVpmxtv3N&branch=master)](https://travis-ci.com/Pas-Kapli/mptp)

## Introduction

The aim of this project is to implement a fast species delimitation method,
based on PTP (Zhang et al. 2013) and the methods and models proposed by
Pons et al. (2006) and Fujisawa et al. (2013). The new tool should:

* have an open source code with an appropriate open source license.
* 64-bit multi-threaded design that handles very large datasets.

We have implemented a tool called mPTP which can handle very large biodiversity
datasets.  It implements a fast method to compute the ML delimitation from an
inferred phylogenetic tree of the samples.  Using MCMC, it also computes the
support values for each clade, which can be used to assess the confidence of
the ML delimitation.

**ML delimitation** mPTP implements two flavours of the point-estimate
solution.  First, it implements the original method from (Zhang et al. 2013)
where all within-species processes are modelled with a single exponential
distribution. mPTP uses a dynamic programming implementation which estimates
the ML delimitation faster and more accurately than the original PTP. The
dynamic programming implementation has similar properties as (Gulek et al.
2010).  See the wiki for more information. The second method assumes a distinct
exponential distribution for the branching events of each of the delimited
species allowing it to fit to a wider range of empirical datasets.

**MCMC method** mPTP generates support values for each clades. They represent
the ratio of the number of samples for which a particular node was in the
between-species process, to the total number of samples. 

## Compilation instructions

**Cloning the repo** Clone the repo and build the executable and the documentation using
the following commands.

```
git clone https://github.com/Pas-Kapli/mptp.git
cd mptp 
./autogen.sh
./configure
make
make install  # as root or sudo make install
```

You will need [GNU Bison](http://www.gnu.org/software/bison/) and
[Flex](http://flex.sourceforge.net/) installed on your system.  Optionally, you
will need the [GNU Scientific Library](http://www.gnu.org/software/gsl/) for
the likelihood ratio test. If it is not available on your system, ratio test
will be disabled.

On a Debian-based Linux system, the three packages can be installed
using the command

`sudo apt-get install libgsl0-dev flex bison`

## Command-line options

General options:

* `--help`
* `--version`
* `--quiet`
* `--tree_show`
* `--multi`
* `--single`
* `--ml`
* `--mcmc INT`
* `--mcmc_sample INT`
* `--mcmc_log`
* `--mcmc_burnin INT`
* `--mcmc_startnull`
* `--mcmc_startrandom`
* `--mcmc_startml`
* `--mcmc_credible REAL`
* `--mcmc_chains INT`
* `--outgroup TAXA`
* `--outgroup_crop`
* `--min_br REAL`
* `--min_br_auto FILENAME`
* `--pvalue REAL`
* `--precision INT`

Input and output options:

* `--tree_file FILENAME`
* `--output_file FILENAME`

Visualization options:

* `--svg_width INT`
* `--svg_fontsize INT`
* `--svg_tipspacing INT`
* `--svg_legend_ratio <0..1>`
* `--svg_nolegend`
* `--svg_marginleft INT`
* `--svg_marginright INT`
* `--svg_margintop INT`
* `--svg_marginbottom INT`
* `--svg_inner_radius INT`

## Usage example

`./mptp --ml --multi --tree_file testTree --output_file out --outgroup A,C --tree_show`

## License and third party licenses

The code is currently licensed under the [GNU Affero General Public License version 3](http://www.gnu.org/licenses/agpl-3.0.en.html).

## Code

    File            | Description
--------------------|----------------
**arch.c**          | Architecture specific code (Mac/Linux).
**auto.c**          | Code for auto-detecting minimum branch length.
**aic.c**           | Code for Bayesian Single- and multi-rate PTP.
**mptp.c**          | Main file handling command-line parameters and executing corresponding parts.
**dp.c**            | Single- and multi-rate DP heuristics for solving the PTP problem.
**fasta.c**         | Code for reading FASTA files.
**lex_rtree.l**     | Lexical analyzer parsing newick rooted trees.
**lex_utree.l**     | Lexical analyzer parsing newick unrooted trees.
**likelihood.c**    | Likelihood rated functions.
**Makefile**        | Makefile.
**maps.c**          | Character mapping arrays for converting sequences to the internal representation.
**output.c**        | Output related files.
**parse_rtree.y**   | Functions for parsing rooted trees in newick format.
**parse_utree.y**   | Functions for parsing unrooted trees in newick format.
**random.c**        | Functions for creating a random delimitation.
**rtree.c**         | Rooted tree manipulation functions.
**svg.c**           | SVG visualization of delimited tree.
**svg_landscape.c** | SVG visualization of likelihood landscape.
**util.c**          | Various common utility functions.
**utree.c**         | Unrooted tree manipulation functions.

## The team

* Paschalia Kapli
* Sarah Lutteropp
* Kassian Kobert
* Pavlos Pavlides
* Jiajie Zhang
* Alexandros Stamatakis
* Tom&aacute;&scaron; Flouri

# References

* Zhang J., Kapli P., Pavlidis P., Stamatakis A. (2013)
**A general species delimitation method with applications to phylogenetic placements.**
*Bioinformatics*, 29(22):2869-2876.
doi:[10.1093/bioinformatics/btt499](http://dx.doi.org/10.1093/bioinformatics/btt499)

* Fujisawa T., Barraclough TG. (2013)
**Delimiting Species Using Single-Locus Data and the Generalized Mixed Yule Coalescent Approach: A Revised Method and Evaluation on Simulated Data Sets.**
*Systematic Biology*, 62(5):707-724.
doi:[10.1093/sysbio/syt033](http://dx.doi.org/10.1093/sysbio/syt033)

* Pons J., Barraclough T., Gomez-Zurita J., Cardoso A., Duran AP, Hazell S., Kamoun S., Sumlin WD, Vogler AP. (2006)
**Sequence-Based Species Delimitation for the DNA Taxonomy of Undescribed Insects.**
*Systematic Biology*, 55(4):595-609.
doi:[10.1080/10635150600852011](http://dx.doi.org/10.1080/10635150600852011)

* Nguyen XV, Epps J., Bailey J. (2010)
**Information Theoretic Measures for Clustering Comparison: Variants, Properties, Normalization and Correction for Chance.**
*Journal of Machine Learning Research*, 11:2837-2854.
[PDF](http://www.jmlr.org/papers/volume11/vinh10a/vinh10a.pdf)

* Gulek M., Toroslu IH. (2010)
**A dynamic programming algorithm for tree-like weighted set packing problem.**
*Information Sciences*, 180(20):3974-3979.
doi:[10.1016/j.ins.2010.06.035](http://dx.doi.org/10.1016/j.ins.2010.06.035)

* Powell JR. (2012)
**Accounting for uncertainty in species delineation during the analysis of environmental DNA sequence data.**
*Methods in Ecology and Evolution*, 3(1):1-11.
doi:[10.1111/j.2041-210X.2011.00122.x](http://dx.doi.org/10.1111/j.2041-210X.2011.00122.x)
