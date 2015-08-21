# Species Delimitation

## Introduction

The aim of this project is to implement a fast species delimitation method,
based on PTP (Zhang et al. 2013) and the methods and models proposed by
Pons et al. (2006) and Fujisawa et al. (2013). The new tool should:

* have an open source code with an appropriate open source license.
* 64-bit multi-threaded design that handles very large datasets.

## Compilation instructions

Currently, the code can be compiled using the included Makefile:

`make`

You will need to have the GNU Scientific Library installed, this can be done
on Ubuntu via `sudo apt-get install libgsl0-dev`

## Command-line options

General options:

* `--help`
* `--version`
* `--tree_show`
* `--ptp_multi`
* `--ptp_single`
* `--outgroup`
* `--quiet`

Input and output options:

* `--tree_file`
* `--output_file`

## Usage example

`./delimit --ptp_multi --tree_file testTree --output_file out --outgroup A,C --tree_show`

## License and third party licenses

The code is currently licensed under the GNU Affero General Public License version 3.

## Code

    File         | Description
-----------------|----------------
**delimit.c**    | Main file handling command-line parameters and executing corresponding parts.
**ptp_multi.c**  | Single- and multi-rate heuristics for solving the PTP assumption.
**Makefile**     | Makefile.
**lex_rtree.l**  | Lexical analyzer parsing newick rooted trees.
**lex_utree.l**  | Lexical analyzer parsing newick unrooted trees.
**util.c**       | Various common utility functions.
**arch.c**       | Architecture specific code (Mac/Linux).
**rtree.c**      | Rooted tree manipulation functions.
**utree.c**      | Unrooted tree manipulation functions.
**parse_rtree.y**| Functions for parsing rooted trees in newick format.
**parse_utree.y**| Functions for parsing unrooted trees in newick format.
**lca_utree.c**  | Naive LCA computation in unrooted trees.
**svg.c**        | SVG visualization.
**score.c**      | Computation of NMI and Kassian's score.

## The team

* Paschalia Kapli
* Sarah Lutteropp
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
