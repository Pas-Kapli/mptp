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

#include "mptp.h"

static double mlsupport_avg(int * mlsupport, double * support, int count)
{
  int i;
  double sum = 0;

  for (i = 0; i < count; ++i)
  {
    if (mlsupport[i] == EVENT_SPECIATION)
      sum += support[i];
    else
      sum += 1 - support[i];
  }

  return sum / count;
}

static void extract_support_recursive(rtree_t * node,
                                      int * index,
                                      double * outbuffer)
{
  if (!node->edge_count) return;

  outbuffer[*index] = node->support;
  *index = *index + 1;

  extract_support_recursive(node->left,  index, outbuffer);
  extract_support_recursive(node->right, index, outbuffer);
}

/* recursively extract support values from a tree into an array */
static int extract_support(rtree_t * root, double * outbuffer)
{
  int index = 0;

  if (!root->edge_count) return -1;

  extract_support_recursive(root, &index, outbuffer);

  return index;
}

static void extract_events_recursive(rtree_t * node,
                                     int * index,
                                     int * outbuffer)
{
  if (!node->edge_count) return;

  outbuffer[*index] = node->event;
  *index = *index + 1;

  extract_events_recursive(node->left,  index, outbuffer);
  extract_events_recursive(node->right, index, outbuffer);
}

static int extract_events(rtree_t * root, int * outbuffer)
{
  int index = 0;

  if (!root->edge_count) return -1;

  extract_events_recursive(root, &index, outbuffer);

  return index;
}

void multichain(rtree_t * root, long method)
{
  long i;
  long * seeds;
  rtree_t * mltree;
  rtree_t ** trees;
  struct drand48_data * rstates;
  double * bayes_min_logl;
  double * bayes_max_logl;

  trees = (rtree_t **)xmalloc((size_t)opt_bayes_chains * sizeof(rtree_t *));
  trees[0] = root;

  /* clone trees in order to have one independent tree per run */
  for (i = 1; i < opt_bayes_chains; ++i)
    trees[i] = rtree_clone(root, NULL);
  mltree = rtree_clone(root,NULL);

  /* allocate memory for storing min and max logl for each chain */
  bayes_min_logl = (double *)xmalloc((size_t)opt_bayes_chains * sizeof(double));
  bayes_max_logl = (double *)xmalloc((size_t)opt_bayes_chains * sizeof(double));

  /* reset to zero */
  memset(bayes_min_logl, 0, (size_t)opt_bayes_chains * sizeof(double));
  memset(bayes_max_logl, 0, (size_t)opt_bayes_chains * sizeof(double));

  /* generate one seed for each run */
  seeds = (long *)xmalloc((size_t)opt_bayes_chains * sizeof(long));
  for (i = 0; i < opt_bayes_chains; ++i)
    seeds[i] = lrand48();
    
  if (opt_bayes_chains == 1)
    seeds[0] = opt_seed;

  /* initialize states for random number generators */
  rstates = (struct drand48_data *)xmalloc((size_t)opt_bayes_chains *
                                           sizeof(struct drand48_data));

  /* initialize a pseudo-random number generator for each chain */
  for (i = 0; i < opt_bayes_chains; ++i)
    srand48_r(seeds[i], rstates+i);

  /* execute each chain sequentially  */
  for (i = 0; i < opt_bayes_chains; ++i)
  {
    dp_init(trees[i]);
    dp_set_pernode_spec_edges(trees[i]);
    if (!opt_quiet)
      fprintf(stdout, "\nBayesian run %ld...\n", i);
    aic_bayes(trees[i],
          method,
          rstates+i,
          seeds[i],
          bayes_min_logl+i,
          bayes_max_logl+i);
    dp_free(trees[i]);

    /* print SVG log-likelihood landscape of current chain given its
       generated seed */
    if (opt_bayes_log)
    {
      svg_landscape(bayes_min_logl[i], bayes_max_logl[i], seeds[i]);
    }

    /* output SVG tree with support values for current chain */
    char * newick = rtree_export_newick(trees[i]);

    if (!opt_quiet)
      fprintf(stdout,
              "Creating tree with support values in %s.%ld.tree ...\n",
              opt_outfile,
              seeds[i]);

    FILE * newick_fp = open_file_ext("tree", seeds[i]);
    fprintf(newick_fp, "%s\n", newick);
    fclose(newick_fp);

    cmd_svg(trees[i], seeds[i]);

    free(newick);
  }

  /* compute the min and max log-l values among all chains */
  double min_logl = bayes_min_logl[0];
  double max_logl = bayes_max_logl[0];
  for (i = 1; i < opt_bayes_chains; ++i)
  {
    if (bayes_min_logl[i] < min_logl) min_logl = bayes_min_logl[i];
    if (bayes_max_logl[i] > max_logl) max_logl = bayes_max_logl[i];
  }

  /* generate the SVG log-likelihood landscape for all chains combined */
  if (!opt_quiet && opt_bayes_log && (opt_bayes_chains > 1))
    fprintf(stdout, "\nPreparing overall log-likelihood landscape ...\n");
  if (opt_bayes_log && (opt_bayes_chains > 1))
    svg_landscape_combined(min_logl, max_logl, opt_bayes_chains, seeds);

  /* free min and max logl arrays */
  free(bayes_min_logl);
  free(bayes_max_logl);

  /* allocate memory for support values */
  double ** support = (double **)xmalloc((size_t)opt_bayes_chains *
                                         sizeof(double *));
  int support_count = 0;
  for (i = 0; i < opt_bayes_chains; ++i)
  {
    support[i] = (double *)xmalloc((size_t)(trees[i]->leaves) * sizeof(double));
    support_count = extract_support(trees[i], support[i]);
    rtree_destroy(trees[i]);
  }

  /* compute ML tree */
  dp_init(mltree);
  dp_set_pernode_spec_edges(mltree);
  dp_ptp(mltree, method);
  int * mlsupport  = (int *)xmalloc((size_t)(mltree->leaves) * sizeof(int));
  if (support_count != extract_events(mltree, mlsupport))
    fatal("Internal error");

  for (i = 0; i < opt_bayes_chains; ++i)
  {
    printf("ML average support based on chain %ld : %.17f\n",
           seeds[i],
           mlsupport_avg(mlsupport, support[i], support_count));
  }

  dp_free(mltree);
  rtree_destroy(mltree);
  free(mlsupport);



  /* compute the standard deviation of each support value given the chains,
     and then compute a consensus average standard deviation for all support
     values */
  double mean, var, stdev, avg_stdev = 0;
  for (i = 0; i < support_count; ++i)
  {
    int j;
    mean = var = stdev = 0;
    for (j = 0; j < opt_bayes_chains; ++j)
      mean += support[j][i];

    mean /= opt_bayes_chains;

    for (j = 0; j < opt_bayes_chains; ++j)
      var += (mean - support[j][i])*(mean - support[j][i]);

    var /= opt_bayes_chains;
    stdev = sqrt(var);

    avg_stdev += stdev;
  }
  avg_stdev /= support_count;

  if (!opt_quiet)
    printf("Average standard deviation of support values among chains: %f\n",
           avg_stdev);

  /* deallocate support values array */
  for (i = 0; i < opt_bayes_chains; ++i)
    free(support[i]);
  free(support);

  /* deallocate all cloned trees (except from the original) */
  free(rstates);
  free(seeds);
  free(trees);
}
