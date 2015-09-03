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

typedef struct cell_s
{
  int coal_edge_count;
  double coal_edgelen_sum;
  double coal_logl;
  double spec_logl;
  char dir;
} cell_t;

static const unsigned char maskup = 1;
static const unsigned char maskdiag = 2;
static const unsigned char maskstop = 4;

static cell_t ** matrix;
static rtree_t ** inner_node_list;

static int all_edge_count = 0;
static double all_edgelen_sum = 0;

static int species_iter = 0;

static double check_lh(cell_t * cell, rtree_t * node, int method)
{
  double coal_logl;
  double spec_logl;

  int edge_count;
  double edgelen_sum;


  if (method == PTP_METHOD_MULTI)
  {
    coal_logl = cell->coal_logl + node->coal_logl;
  }
  else
  {
    /* compute coalescent edges */
    edge_count = cell->coal_edge_count + node->edge_count;
    edgelen_sum = cell->coal_edgelen_sum + node->edgelen_sum;
    assert(edge_count >= 0);
    coal_logl = loglikelihood(edge_count, edgelen_sum);
  }

  /* compute speciation edge count and edge length sum */
  edge_count = all_edge_count - cell->coal_edge_count - node->edge_count;
  edgelen_sum = all_edgelen_sum - cell->coal_edgelen_sum - node->edgelen_sum;

  spec_logl = loglikelihood(edge_count, edgelen_sum);

  return (spec_logl + coal_logl);
}

static void update_cell(cell_t * from, cell_t * to, rtree_t * node, int method)
{
  int edge_count;
  double edgelen_sum;

  if (method == PTP_METHOD_MULTI)
  {
    to->coal_logl = from->coal_logl + node->coal_logl; 
  }
  else
  {
    edge_count = from->coal_edge_count + node->edge_count;
    edgelen_sum = from->coal_edgelen_sum + node->edgelen_sum;
    to->coal_logl = loglikelihood(edge_count, edgelen_sum);
  }

  /* compute speciation edge count and edge length sum */
  edge_count = all_edge_count - from->coal_edge_count - node->edge_count;
  edgelen_sum = all_edgelen_sum - from->coal_edgelen_sum - node->edgelen_sum;

  to->spec_logl = loglikelihood(edge_count, edgelen_sum);


  to->coal_edge_count = all_edge_count - edge_count;
  to->coal_edgelen_sum = all_edgelen_sum - edgelen_sum;

}

static double cell_logl(cell_t * cell)
{
  return (cell->coal_logl + cell->spec_logl);
}

static void dp_knapsack_recurse(rtree_t * node, int method)
{
  static int i = 1;

  int j,c;

  if (!node->left) return;

  /* bottom-up traversal */
  dp_knapsack_recurse(node->left,  method);
  dp_knapsack_recurse(node->right, method);

  /* allocate current row of matrix */
  matrix[i] = (cell_t *)xmalloc(((size_t)i+1) * sizeof(cell_t));

  /* traverse all inner nodes of subtree except self and copy cell from above */
  for (j = 0; j <  node->leaves-1; ++j)
  {
    memcpy(&matrix[i][j], &matrix[i-1][j], sizeof(cell_t));
    matrix[i][j].dir = maskup;
  }

  /* iterate through the rest of cells */
  for (j = node->leaves-1; j < i; ++j)
  {
    c = node->leaves - 1;
    cell_t * diag = &matrix[i-c][j-c];
    cell_t * top  = &matrix[i-1][j];

    if (check_lh(diag, node, method) > (top->coal_logl + top->spec_logl))
    {
      update_cell(diag, &matrix[i][j], node, method);
      matrix[i][j].dir = maskdiag;
    }
    else
    {
      memcpy(&matrix[i][j], top, sizeof(cell_t));
      matrix[i][j].dir = maskup;
    }
  }

  /* last cell is this */
  c = node->leaves-1;
  update_cell(&matrix[i-c][j-c], &matrix[i][i], node, method);
  matrix[i][i].dir = maskdiag;

  inner_node_list[i] = node;

  ++i;
}

static void ks_init(rtree_t * tree)
{
  if (tree->left)  ks_init(tree->left);
  if (tree->right) ks_init(tree->right);

  tree->coal_logl = loglikelihood(tree->edge_count,
                                  tree->edgelen_sum);
}

static int ks_backtrack(int i, int j)
{
  int species_count = inner_node_list[i]->leaves;

  while (matrix[i][j].dir != maskstop)
  {
    if (matrix[i][j].dir == maskup)
    {
      inner_node_list[i]->event = EVENT_SPECIATION;
      --i;
    }
    else
    {
      species_count = species_count - inner_node_list[i]->leaves + 1;
      inner_node_list[i]->event = EVENT_COALESCENT;
      int c = inner_node_list[i]->leaves - 1;
      i -= c;
      j -= c;
      assert(i >= 0 && j >= 0);
    }
  }

  return species_count;
}

static void ks_delimit(rtree_t * node, FILE * out, int * warning)
{
  if (node->event == EVENT_COALESCENT)
  {
    species_iter++;

    fprintf(out, "\nSpecies %d:\n", species_iter);
    rtree_print_tips(node,out);
  }
  else
  {
    if (node->length <= opt_minbr && node->parent) *warning = 1;
    ks_delimit(node->left,out,warning);
    ks_delimit(node->right,out,warning);
  }
}

void dp_knapsack(rtree_t * root, int method)
{
  int inner_nodes_count = root->leaves - 1;
  int i;
  double score;

  ks_init(root);

  /* allocate matrix */
  matrix = (cell_t **)xmalloc(((size_t)inner_nodes_count+1) * sizeof(cell_t *));
  matrix[0] = (cell_t *)xmalloc(sizeof(cell_t));
  inner_node_list = (rtree_t **)xmalloc(((size_t)inner_nodes_count+1) * 
                                        sizeof(rtree_t *));
  inner_node_list[0] = NULL;

  /* initialize global variables */
  all_edge_count = root->edge_count;
  all_edgelen_sum = root->edgelen_sum;

  /* initialize */
  matrix[0][0].coal_logl = 0;
  matrix[0][0].spec_logl = loglikelihood(all_edge_count, all_edgelen_sum);
  matrix[0][0].coal_edge_count = 0;
  matrix[0][0].coal_edgelen_sum = 0;
  matrix[0][0].dir = maskstop;


  /* bottom-up traversal and computation */
  dp_knapsack_recurse(root, method);

  int best_index = 0;
  score = cell_logl(&matrix[inner_nodes_count][0]);
  for (i = 0; i <= inner_nodes_count; ++i)
  {
    if (cell_logl(&matrix[inner_nodes_count][i]) > score)
    {
      score = cell_logl(&matrix[inner_nodes_count][i]);
      best_index = i;
    }
  }

  if (!opt_quiet)
  {
    fprintf(stdout, 
           "Number of edges greater than minimum branch length: %d / %d\n", 
           root->edge_count,
           2 * root->leaves - 2);
    fprintf(stdout,"Null-model score: %.6f\n", root->coal_logl);
    if (method == PTP_METHOD_SINGLE)
      fprintf(stdout, "Best score for single coalescent rate: %.6f\n", score); 
    else
      fprintf(stdout, "Best score for multi coalescent rate: %.6f\n", score);
  }

  /* backtrack and mark nodes with coalescent/speciation events */
  int species_count = ks_backtrack(inner_nodes_count, best_index);
  assert(species_count > 0);

  /* do a Likelihood Ratio Test (lrt) and return the computed p-value */
  double pvalue = -1;
  int lrt_pass;
  lrt_pass = lrt(root->coal_logl,
                 score,
                 (method == PTP_METHOD_MULTI) ?
                   (unsigned int)species_count : 1,
                 &pvalue);

  if (!opt_quiet)
    fprintf(stdout,"LRT computed p-value: %.6f\n", pvalue);

  FILE * out = open_file_ext(".txt");

  if (!opt_quiet)
    fprintf(stdout, "Writing delimitation file %s.txt ...\n", opt_outfile);

  /* write information about delimitation to file */
  output_info(out,
              method,
              root->coal_logl,
              score,
              pvalue,
              lrt_pass,
              root,
              lrt_pass ?
                species_count : 1);

  /* if LRT passed, then back-track the DP table and print the delimitation,
     otherwise print the null-model (one single species) */
  if (lrt_pass)
  {
    int warning = 0;
    ks_delimit(root,out, &warning);
    if (warning)
      fprintf(stderr,"WARNING: A speciation edge is smaller than the specified "
                     "minimum branch length.\n");
  }
  else
  {
    species_iter = 1;
    if (!opt_quiet)
      fprintf(stdout, "LRT failed -- null-model is preferred and printed\n");
    fprintf(out,"Species 1:\n");
    rtree_print_tips(root,out);
  }

  fclose(out);

  if (!opt_quiet)
    fprintf(stdout, "Number of delimited species: %d\n", species_iter);

  if (root->edge_count == 0)
    fprintf(stderr, "WARNING: The tree has no edges > %f. "
                    "All edges have been ignored. \n", opt_minbr);

  /* deallocate DP matrix */
  for (i = 0; i <= inner_nodes_count; ++i)
    free(matrix[i]);
  free(inner_node_list);
  free(matrix);
}
