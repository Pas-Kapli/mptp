/*
    Copyright (C) 2015 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "delimit.h"

void multichain(rtree_t * root, int method, prior_t * prior)
{
  long i;
  long * seeds;
  rtree_t ** trees;
  struct drand48_data * rstates;
  double * bayes_min_logl;
  double * bayes_max_logl;

  trees = (rtree_t **)xmalloc(opt_bayes_chains * sizeof(rtree_t *));
  trees[0] = root;

  /* clone trees in order to have one independent tree per run */
  for (i = 1; i < opt_bayes_chains; ++i)
    trees[i] = rtree_clone(root, NULL);

  bayes_min_logl = (double *)xmalloc(opt_bayes_chains * sizeof(double));
  bayes_max_logl = (double *)xmalloc(opt_bayes_chains * sizeof(double));

  memset(bayes_min_logl, 0, opt_bayes_chains * sizeof(double));
  memset(bayes_max_logl, 0, opt_bayes_chains * sizeof(double));

  /* generate a seed for each run */
  seeds = (long *)xmalloc(opt_bayes_chains * sizeof(long));
  for (i = 0; i < opt_bayes_chains; ++i)
    seeds[i] = lrand48();

  /* initialize states for random number generators */
  rstates = (struct drand48_data *)xmalloc(opt_bayes_chains *
                                           sizeof(struct drand48_data));
  for (i = 0; i < opt_bayes_chains; ++i)
    srand48_r(seeds[i], rstates+i);

  for (i = 0; i < opt_bayes_chains; ++i)
  {
    assert(trees[i]->parent == NULL);
    dp_init(trees[i]);
    dp_set_pernode_spec_edges(trees[i]);
    bayes(trees[i], method, prior, rstates+i, seeds[i], bayes_min_logl+i, bayes_max_logl+i);
    dp_free(trees[i]);
    svg_landscape(bayes_min_logl[i], bayes_max_logl[i], seeds[i]);
    
    /* output tree with support values */
    char * newick = rtree_export_newick(trees[i]);

    if (!opt_quiet)
      fprintf(stdout,
              "Creating tree with support values in %s.%ld.tree ...\n",
              opt_outfile,
              seeds[i]);

    FILE * newick_fp = open_file_ext("tree", seeds[i]);
    fprintf(newick_fp, "%s\n", newick);
    fclose(newick_fp);

    cmd_svg(trees[i]);

    free(newick);
  }

  double min_logl = bayes_min_logl[0];
  double max_logl = bayes_min_logl[0];
  for (i = 1; i < opt_bayes_chains; ++i)
  {
    if (bayes_min_logl[i] < min_logl) min_logl = bayes_min_logl[i];
    if (bayes_max_logl[i] > max_logl) max_logl = bayes_max_logl[i];
  }

  if (!opt_quiet)
    fprintf(stdout, "\nPreparing overall log-likelihood landscape ...\n");
  svg_landscape_combined(min_logl, max_logl, opt_bayes_chains, seeds);

  free(bayes_min_logl);
  free(bayes_max_logl);

  /* deallocate all cloned trees (except from the original) */
  for (i = 0; i < opt_bayes_chains; ++i)
    rtree_destroy(trees[i]);
  free(rstates);
  free(seeds);
  free(trees);
}
