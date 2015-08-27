/*
    Copyright (C) 2015 Tomas Flouri, Sarah Lutteropp

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

static unsigned int species_iter = 0;

static double loglikelihood(int edge_count, double edgelen_sum)
{
  assert(edge_count >= 0);

  /* TODO: Find out whether edgelen_sum == 0 makes sense */
  if (edge_count == 0 || edgelen_sum == 0) return 0;
  
  return edge_count * (log(edge_count) - 1 - log(edgelen_sum));
}

static void dp_recurse(rtree_t * node, int method, prior_t * prior)
{
  int k,j;

  /* bottom-up recursion */

  if (node->left)  dp_recurse(node->left,  method, prior);
  if (node->right) dp_recurse(node->right, method, prior);

  /*                u_vec
                *
               / \
              /   \
     v_vec   *     *  w_vec    */
  
  dp_vector_t * u_vec = node->vector;

  double spec_logl = loglikelihood(node->spec_edge_count,
                                   node->spec_edgelen_sum) +
                     prior_score(1,prior);

  u_vec[0].spec_edgelen_sum = 0;
  u_vec[0].score_multi = node->coal_logl + spec_logl;
  u_vec[0].score_single = node->coal_logl + spec_logl;
  u_vec[0].coal_multi_logl = node->coal_logl;
  u_vec[0].species_count = 1;
  u_vec[0].filled = 1;

  if (!node->left) return;

  dp_vector_t * v_vec = node->left->vector;
  dp_vector_t * w_vec = node->right->vector;

  assert(node->spec_edge_count >= 0);

  int u_edge_count = 0;
  double u_edgelen_sum = 0;

  /* check whether edges (u,v) and (u,w) are > min branch length */
  if (node->left->length > opt_minbr)
  {
    u_edge_count++;
    u_edgelen_sum += node->left->length;
  }
  if (node->right->length > opt_minbr)
  {
    u_edge_count++;
    u_edgelen_sum += node->right->length;
  }

  for (j = 0; j <= node->left->edge_count; ++j)
  {
    for (k = 0; k <= node->right->edge_count; ++k)
    {
      /* if at least one of the two entries is not valid/filled, skip */
      if (!v_vec[j].filled || !w_vec[k].filled) continue;

      int i = j + k + u_edge_count;

      /* set the number of species - needed for prior information */
      unsigned int species_count = v_vec[j].species_count +
                                   w_vec[k].species_count;

      /* compute multi-rate coalescent log-likelihood */
      double coal_multi_logl = v_vec[j].coal_multi_logl +
                               w_vec[k].coal_multi_logl;
      
      /* compute coalescent edge count and length sum of subtree u */
      double u_spec_edgelen_sum = v_vec[j].spec_edgelen_sum +
                                  w_vec[k].spec_edgelen_sum +
                                  u_edgelen_sum;
      int coal_edge_count = node->edge_count - i;            /* change to int */
      double coal_edgelen_sum = node->edgelen_sum - u_spec_edgelen_sum;
                                                  
                                                  
      /* compute single-rate coalescent log-likelihood */
      double coal_single_logl = loglikelihood(coal_edge_count,coal_edgelen_sum);

      /* compute total speciation log-likelihood */
      double spec_edgelen_sum = node->spec_edgelen_sum +
                                u_edgelen_sum +
                                v_vec[j].spec_edgelen_sum +
                                w_vec[k].spec_edgelen_sum;

      int spec_edge_count  = node->spec_edge_count + i;
      assert(species_count > 0);
      spec_logl = loglikelihood(spec_edge_count,spec_edgelen_sum) +
                  prior_score(species_count,prior);

      
      /* compute single- and multi-rate scores */
      double score_multi = coal_multi_logl + spec_logl;
      double score_single = coal_single_logl + spec_logl;
      double score = score_multi;
      double best_score = u_vec[i].score_multi;

      if (method == PTP_METHOD_SINGLE)
      {
        score = score_single;
        best_score = u_vec[i].score_single;
      }

      if (!u_vec[i].filled || score > best_score)
      {
        u_vec[i].score_multi = score_multi;
        u_vec[i].score_single = score_single;
        u_vec[i].spec_edgelen_sum = u_spec_edgelen_sum;
        u_vec[i].coal_multi_logl = coal_multi_logl;
        u_vec[i].vec_left = j;
        u_vec[i].vec_right = k;
        u_vec[i].species_count = species_count;
        u_vec[i].filled = 1;
      }

    }
  }
}

static void print_subtree_leaves(rtree_t * root)
{
  if (root->left)
    print_subtree_leaves(root->left);
  if (root->right)
    print_subtree_leaves(root->right);

  if (!root->left && !root->right)
    fprintf(stdout, "%s\n", root->label);
}

static void backtrack(rtree_t * node,
                      int index,
                      bool *warning_minbr)

{
  dp_vector_t * vec = node->vector;

  if ((vec[index].vec_left != -1) && (vec[index].vec_right != -1))
  {
    node->event = EVENT_SPECIATION;

    if (node->length <= opt_minbr && node->parent) *warning_minbr = true;

    backtrack(node->left, vec[index].vec_left, warning_minbr);
    backtrack(node->right,vec[index].vec_right,warning_minbr);
  }
  else
  {
    species_iter++;
    node->event = EVENT_COALESCENT;

    if (!opt_quiet)
    {
      fprintf(stdout, "\nSpecies %d:\n", species_iter);
      print_subtree_leaves(node);
    }
  }
}

static void print_null_model(rtree_t * tree)
{
  printf("Species 1:\n");
  print_subtree_leaves(tree);
}

static int lrt(double nullmodel_logl,
               double ptp_logl,
               unsigned int df,
               double * pvalue)
{
  double diff = 2*(ptp_logl - nullmodel_logl);

  /* http://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chdtr.html */
  *pvalue = 1 - gsl_cdf_chisq_P(diff,df);

  if ((*pvalue) > opt_pvalue)
    return 0;
  
  return 1;
}

void dp_ptp(rtree_t * tree, int method, prior_t * prior)
{
  int i;

  /* reset species counter */
  species_iter = 0;

  /* fill DP table */
  dp_recurse(tree, method, prior);

  double max = 0;
  int best_index = 0;

  /* obtain best entry in the root DP table */
  dp_vector_t * vec = tree->vector;
  if (method == PTP_METHOD_MULTI)
    for (i = 0; i < tree->edge_count; i++)
    {
      if (max < vec[i].score_multi)
      {
        max = vec[i].score_multi;
        best_index = i;
      }
    }
  else
    for (i = 0; i < tree->edge_count; i++)
    {
      if (max < vec[i].score_single)
      {
        max = vec[i].score_single;
        best_index = i;
      }
    }

  /* output some statistics */
  if (!opt_quiet)
  {
    printf("Num big enough edges in tree: %d\n",
      tree->edge_count);
    printf("Score Null Model: %.6f\n", tree->coal_logl);
    printf("Best score found single: %.6f\n", vec[best_index].score_single);
    printf("Best score found multi: %.6f\n", vec[best_index].score_multi);
  }

  /* do a Likelihood Ratio Test (lrt) and return the computed p-value */
  double pvalue = -1;
  int lrt_pass;
  if (method == PTP_METHOD_MULTI)
  {
    lrt_pass = lrt(tree->coal_logl, 
                   vec[best_index].score_multi,
                   vec[best_index].species_count,
                   &pvalue);
  }
  else
  {
    lrt_pass = lrt(tree->coal_logl,
                   vec[best_index].score_single,
                   1,
                   &pvalue);
  }

  if (!opt_quiet)
    printf("Computed P-value: %.6f\n", pvalue);

  /* if LRT passed, then back-track the DP table and print the delimitation,
     otherwise print the null-model (one single species) */
  if (lrt_pass)
  {
    bool warning_minbr = false;
    backtrack(tree, best_index, &warning_minbr);
    if (!opt_quiet && warning_minbr)
      printf("WARNING: A speciation edge in the result is too small.\n");
  }
  else
  {
    species_iter = 1;
    if (!opt_quiet)
    {
      printf("The Null Model is the preferred one.\n");
      print_null_model(tree);
    }
  }
  if (!opt_quiet)
  {
    printf("Number of delimited species: %d\n", species_iter);
  }
  if (tree->edge_count == 0 && !opt_quiet)
  {
    fprintf(stderr, "WARNING: The tree has no edges > %f. All edges have been ignored. \n", opt_minbr);
  }
}

void dp_init(rtree_t * tree)
{
  int i;

  if (tree->left)  dp_init(tree->left);
  if (tree->right) dp_init(tree->right);

  // TODO: Check whether this is the best way to handle those
  //   nasty zero-length edges.

  tree->vector = calloc((size_t)(tree->edge_count + 1), sizeof(dp_vector_t));

  for (i = 0; i <= tree->edge_count; i++)
  {
    tree->vector[i].vec_left  = -1;
    tree->vector[i].vec_right = -1;
  }

  assert(tree->edge_count >= 0);

  tree->coal_logl = loglikelihood(tree->edge_count,
                                  tree->edgelen_sum);
}

void dp_free(rtree_t * tree)
{
  if (tree->left)  dp_free(tree->left);
  if (tree->right) dp_free(tree->right);

  if (tree->vector) free(tree->vector);
}

void dp_set_pernode_spec_edges(rtree_t * node)
{
  if (!node) return;

  node->spec_edge_count = 0;
  node->spec_edgelen_sum = 0;

  /* for each node set spec_edge_count (and spec_edgelen_sum) as the count
     (or sum) of edges (edge-lengths) of all direct child edges of 
     nodes on the path to root excluding the current node */
  if (node->parent)
  {
    node->spec_edge_count = node->parent->spec_edge_count;
    node->spec_edgelen_sum = node->parent->spec_edgelen_sum;

    double len = node->parent->left->length;
    if (len > opt_minbr) 
    {
      node->spec_edge_count++; 
      node->spec_edgelen_sum += len;
    }

    len = node->parent->right->length;
    if (len > opt_minbr) 
    {
      node->spec_edge_count++; 
      node->spec_edgelen_sum += len;
    }
  }

  dp_set_pernode_spec_edges(node->left);
  dp_set_pernode_spec_edges(node->right);
}

