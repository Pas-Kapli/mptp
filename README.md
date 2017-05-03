# Species Delimitation

[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)
[![Build Status](https://travis-ci.org/Pas-Kapli/mptp.svg?branch=master)](https://travis-ci.com/Pas-Kapli/mptp)

## Introduction

The aim of this project is to implement a fast species delimitation method,
based on PTP (Zhang et al. 2013). The new tool should:

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
2010).  See the [wiki](https://github.com/Pas-Kapli/mptp/wiki) for more
information. The second method assumes a distinct exponential distribution for
the branching events of each of the delimited species allowing it to fit to a
wider range of empirical datasets.

**MCMC method** mPTP generates support values for each clades. They represent
the ratio of the number of samples for which a particular node was in the
between-species process, to the total number of samples. 

## Compilation instructions

**Cloning the repo** Clone the repo and build the executable and the documentation using
the following commands.

```bash
git clone https://github.com/Pas-Kapli/mptp.git
cd mptp 
./autogen.sh
./configure
make
make install  # as root, or run sudo make install
```

You will need [GNU Bison](http://www.gnu.org/software/bison/) and
[Flex](http://flex.sourceforge.net/) installed on your system.  When using the
cloned repository version, you will also need
[autoconf](https://www.gnu.org/software/autoconf/autoconf.html) and
[automake](https://www.gnu.org/software/automake/) installed. Optionally, you
will need the [GNU Scientific Library](http://www.gnu.org/software/gsl/) for
the likelihood ratio test. If it is not available on your system, ratio test
will be disabled.

On a Debian-based Linux system, the four packages can be installed
using the command

```bash
sudo apt-get install libgsl0-dev flex bison autotools-dev autoconf
```

Optionally, you can install the bash auto-completion for mptp. To do that,
replace the `./configure` step above with
```bash
./configure --with-bash-completions=DIR
```
where `DIR` is the directory where bash autocompletion is stored. You can use
`pkg-config` as follows:
```bash
./configure --with-bash-completions=`pkg-config --variable=completionsdir bash-completion`
```

**Source distribution** To download the source distribution from a
[release](https://github.com/Pas-Kapli/mptp/releases) and build the executable
and the documentation, use the following commands:

```bash
wget https://github.com/Pas-Kapli/mptp/releases/download/v0.2.2/mptp-src-0.2.2.tar.gz
tar zxvf mptp-src-0.2.2.tar.gz
cd mptp-src-0.2.2
./configure
make
make install  # as root, or run sudo make install
```

Note that, similarly to cloning the repository, you will need [GNU
Bison](http://www.gnu.org/software/bison/) and
[Flex](http://flex.sourceforge.net/) installed on your system, and optionally,
the [GNU Scientific Library](http://www.gnu.org/software/gsl/).  However, you
do not need [autoconf](https://www.gnu.org/software/autoconf/autoconf.html) and
[automake](https://www.gnu.org/software/automake/) installed (note the missing `./autogen`).
See also the notes for installing the bash auto-completition, as described in
the *Cloning the repo* section.


**Binary distribution** Starting with version 0.2.0, binary distribution files
(.tar.gz) for GNU/Linux on x86-64 containing pre-compiled binaries as well as
the documentation (man and pdf files) will be made available as part of each
[release](https://github.com/Pas-Kapli/mptp/releases). The included executables
currently are not compiled with [`libgsl`](http://www.gnu.org/software/gsl/)
support. This means, Likelihood Ratio Test (LRT) is disabled for the
single-rate PTP model. However, we intend to implement dynamic loading for
`libgsl` and therefore this issue will disappear in the next releases. Until then, please
consider compiling from source in order to enable `libgsl`.

To use the pre-compiled binary, download the appropriate executable for your
system using the following commands if you are using a Linux system:

```bash
wget https://github.com/Pas-Kapli/mptp/releases/download/v0.2.2/mptp-0.2.2-linux-x86_64.tar.gz
tar zxvf mptp-0.2.2-linux-x86_64.tar.gz
```

You will now have the binary distribution in a folder called
`mptp-0.2.2-linux-x86_64` in which you will find three subfolders `bin`, `man`
and `doc`. We recommend making a copy or a symbolic link to the mptp binary
`bin/mptp` in a folder included in your `$PATH`, and a copy or a symbolic link
to the mptp man page `man/mptp.1` in a folder included in your `$MANPATH`. The
PDF version of the manual is available in `doc/mptp_manual.pdf`.



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
* `--mcmc_runs INT`
* `--outgroup TAXA`
* `--outgroup_crop`
* `--minbr REAL`
* `--minbr_auto FILENAME`
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

```bash
mptp --ml --multi --tree_file testTree --output_file out --outgroup A,C --tree_show
mptp --mcmc 50000000 --multi --mcmc_sample 1000000 --mcmc_burnin 1000000 --tree_file tree.newick --output_file out
```

## Documentation

If `mptp` was installed according to the [Compilation
instructions](https://github.com/Pas-Kapli/mptp#compilation-instructions) you
can access the man pages by:

```bash
man mptp
```

A comprehensive documentation is also available in the [wiki](https://github.com/Pas-Kapli/mptp/wiki).

## License and third party licenses

The code is currently licensed under the [GNU Affero General Public License version 3](http://www.gnu.org/licenses/agpl-3.0.en.html).

## Code

    File            | Description
--------------------|----------------
**arch.c**          | Architecture specific code (Mac/Linux).
**auto.c**          | Code for auto-detecting minimum branch length.
**aic.c**           | Code for Bayesian Single- and multi-rate PTP.
**mptp.c**          | Main file handling command-line parameters and executing corresponding parts.
**mptp.h**          | MPTP Header file.
**dp.c**            | Single- and multi-rate DP heuristics for solving the PTP problem.
**fasta.c**         | Code for reading FASTA files.
**lex_rtree.l**     | Lexical analyzer parsing newick rooted trees.
**lex_utree.l**     | Lexical analyzer parsing newick unrooted trees.
**likelihood.c**    | Likelihood rated functions.
**Makefile.am**     | Automake file for generating Makefile.in.
**maps.c**          | Character mapping arrays for converting sequences to the internal representation.
**multirun.c**      | Functions to execute multiple MCMC runs and compute ASD of support values.
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
