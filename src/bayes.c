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

#include "mptp.h"

typedef struct density_s
{
  double logl;
  long species_count;
} density_t;


static rtree_t ** crnodes;
static rtree_t ** snodes;

static int crnodes_count = 0;
static int snodes_count = 0;

static long accept_count = 0;
static FILE * fp_log = NULL;

static long species_count = 0;

static density_t * densities = NULL;

static void mcmc_log(double logl, long sc)
{
  if (opt_bayes_log)
    fprintf(fp_log, "%f,%ld\n", logl, sc);
}

int cb_desc(const void * va, const void * vb)
{
  const density_t * a = va;
  const density_t * b = vb;

  if (a->logl - b->logl < 0)
    return 1;
  else if (a->logl - b->logl > 0)
    return -1;

  return 0;
}

static void bayes_init(rtree_t * root, long seed)
{
  long i;

  crnodes = (rtree_t **)xmalloc(root->leaves*sizeof(rtree_t *));
  snodes = (rtree_t **)xmalloc(root->leaves*sizeof(rtree_t *));

  crnodes_count = 0;
  snodes_count = 0;
  accept_count = 0;

  densities = (density_t *)xmalloc((root->leaves+1)*sizeof(density_t));
  memset(densities, 0, (root->leaves+1) * sizeof(density_t));
  for (i = 0; i < root->leaves+1; ++i)
    densities[i].species_count = i;

  /* open log file */
  if (opt_bayes_log)
    fp_log = open_file_ext("log", seed);
}

static void init_null(rtree_t * root)
{
  int i;

  rtree_t ** inner_node_list = (rtree_t **)xmalloc((root->leaves-1) *
                                                   sizeof(rtree_t *));
  rtree_query_innernodes(root, inner_node_list);

  /* start bayes analysis from null model */
  for (i = 0; i < root->leaves - 1; ++i)
    inner_node_list[i]->event = EVENT_COALESCENT;
  free(inner_node_list);
}

static void bayes_stats_init(rtree_t * root)
{
  int i;

  rtree_t ** inner_node_list = (rtree_t **)xmalloc((root->leaves-1) *
                                                   sizeof(rtree_t *));
  rtree_query_innernodes(root, inner_node_list);

  for (i = 0; i < root->leaves - 1; ++i)
  {
    if (inner_node_list[i]->event == EVENT_COALESCENT)
      inner_node_list[i]->speciation_start = -1;
    else
      inner_node_list[i]->speciation_start = opt_bayes_burnin-1;

    inner_node_list[i]->speciation_count = 0;
  }

  free(inner_node_list);
}

static void hpd(long n, FILE * fp)
{
  long i;
  long min, max;
  double densities_sum = 0;
  double acc_sum = 0;
  long * indices = NULL;

  indices = (long *)xmalloc((n+2)*sizeof(long));
  memset(indices, 0, (n+2) * sizeof(long));

  for (i = 1; i <= n; ++i)
    densities_sum += densities[i].logl;

  max = 0; min = n+1;
  for (i = 1; i <= n; ++i)
  {
    acc_sum += densities[i].logl;
    indices[densities[i].species_count] = 1;

    if (densities[i].species_count < min)
      min = densities[i].species_count;

    if (densities[i].species_count > max)
      max = densities[i].species_count;

    if (acc_sum / densities_sum >= opt_bayes_credible)
      break;
  }

  fprintf(fp, "CCI (%ld,%ld)\n", min, max);
  if (!opt_quiet)
    fprintf(stdout, "CCI (%ld,%ld)\n", min, max);

  
  fprintf(fp, "HPD ");
  if (!opt_quiet)
    printf("HPD ");
  for (i = 1; i <= n+1; ++i)
  {
    if (indices[i] == 1 && indices[i-1] == 0)
    {
      fprintf(fp, "(%ld,", i);
      if (!opt_quiet)
        printf("(%ld,", i);
    }
    if (indices[i] == 0 && indices[i-1] == 1)
    {
      fprintf(fp, "%ld) ", i-1);
      if (!opt_quiet)
        printf("%ld) ", i-1);
    }
  }
  fprintf(fp,"\n");
  if (!opt_quiet)
    printf("\n");
  free(indices);

}

static void bayes_finalize(rtree_t * root,
                           double bayes_min_logl,
                           double bayes_max_logl,
                           long seed)
{
  long i;

  if (!opt_quiet)
    printf ("Maximum log-likelihood observed in bayesian run: %f\n", bayes_max_logl);

  /* write support values to all nodes */
  rtree_t ** inner_node_list = (rtree_t **)xmalloc((root->leaves-1) *
                                                   sizeof(rtree_t *));
  rtree_query_innernodes(root, inner_node_list);

  for (i = 0; i < root->leaves - 1; ++i)
  {
    if (inner_node_list[i]->speciation_start != -1)
      inner_node_list[i]->speciation_count = inner_node_list[i]->speciation_count +
                                             opt_bayes_runs - 
                                             inner_node_list[i]->speciation_start;

    inner_node_list[i]->support = inner_node_list[i]->speciation_count / 
                                  (double)(opt_bayes_runs-opt_bayes_burnin+1);
  }

  free(inner_node_list);
  free(crnodes);
  free(snodes);

  if (opt_bayes_log)
  {
    if (!opt_quiet)
      fprintf(stdout, "Log written in %s.%ld.log ...\n", opt_outfile, seed);

    fclose(fp_log);
  }

  FILE * fp_stats = open_file_ext("stats", seed);

  double densities_sum = 0;
  for (i = 1; i <= root->leaves; ++i)
    densities_sum += densities[i].logl;
    
  for (i = 1; i <= root->leaves; ++i)
  {
    fprintf(fp_stats,
            "%ld,%f\n",
            i,
            (densities[i].logl/densities_sum)*100);
  }

  /* compute a HPD */
  qsort(densities+1, root->leaves, sizeof(density_t), cb_desc);
  hpd(root->leaves, fp_stats);

  if (!opt_quiet)
    fprintf(stdout,
            "Statistics written in %s.%ld.stats ...\n",
            opt_outfile,
            seed);

  fclose(fp_stats);
  free(densities);
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

static void backtrack_random(rtree_t * node,
                             bool *warning_minbr)

{

  node->bayes_slot = -1;

  if (node->event == EVENT_SPECIATION)
  {
    if (node->length <= opt_minbr && node->parent) *warning_minbr = true;

    backtrack_random(node->left,  warning_minbr);
    backtrack_random(node->right, warning_minbr);

    /* add to list of speciation nodes only if its two direct descendents
       are coalescent roots and also the subtree at node has at least one
       branch length greater than minbr */
    if ((node->left->event == EVENT_COALESCENT) && 
        (node->right->event == EVENT_COALESCENT) &&
        (node->edge_count))
    {
      node->bayes_slot = snodes_count;
      snodes[snodes_count++] = node;
    }

  }
  else
  {
    
    node->event = EVENT_COALESCENT;

    /* add to list of coalescent roots in case it is not a tip AND if
       the subtree rooted at node has at least one edge longer than minbr */
    if (node->edge_count)
    {
      node->bayes_slot = crnodes_count;
      crnodes[crnodes_count++] = node;
    }
  }
}

static void backtrack(rtree_t * node,
                      int index,
                      bool *warning_minbr)

{
  dp_vector_t * vec = node->vector;

  node->bayes_slot = -1;

  if ((vec[index].vec_left != -1) && (vec[index].vec_right != -1))
  {
    node->event = EVENT_SPECIATION;

    if (node->length <= opt_minbr && node->parent) *warning_minbr = true;

    backtrack(node->left, vec[index].vec_left, warning_minbr);
    backtrack(node->right,vec[index].vec_right,warning_minbr);

    /* add to list of speciation nodes only if its two direct descendents
       are coalescent roots and also the subtree at node has at least one
       branch length greater than minbr */
    if ((node->left->event == EVENT_COALESCENT) && 
        (node->right->event == EVENT_COALESCENT) &&
        (node->edge_count))
    {
      node->bayes_slot = snodes_count;
      snodes[snodes_count++] = node;
    }

  }
  else
  {
    node->event = EVENT_COALESCENT;

    /* add to list of coalescent roots in case it is not a tip AND if
       the subtree rooted at node has at least one edge longer than minbr */
    if (node->edge_count)
    {
      node->bayes_slot = crnodes_count;
      crnodes[crnodes_count++] = node;
    }
  }
}

static void speciate(unsigned int r)
{
  /*            CR                         S
                *                          *
               / \            ->          / \
              /   \                      /   \
          C  *     *  C             CR  *     *  CR            */


  /* select the coalescent root at position r and split it into 
     two coalescent root nodes */

  rtree_t * node = crnodes[r];
  
  /* move the last node of the list to the position of the node
     we just used */
  if (r != (unsigned int)(crnodes_count-1))
  {
    crnodes[r] = crnodes[crnodes_count-1];
    crnodes[r]->bayes_slot = r;
  }
  --crnodes_count;

  /* eliminate parent from snodes if both its children were coalescent
     roots, i.e. we had the case below:

              S                               S
              *                               *
             / \                             / \
            /   \                           /   \
       CR  *     *  CR         ->      CR  *     *  S
                / \                             / \
               /   \                           /   \
           C  *     *  C                  CR  *     *  CR
 
  */
  if (node->parent && 
      node->parent->left->event == EVENT_COALESCENT && 
      node->parent->right->event == EVENT_COALESCENT)
  {
    assert(node->parent->bayes_slot != -1);
    assert(node->edge_count);

    /* perform the following only if the parent is not the last node
       in the list */
    if (node->parent->bayes_slot != snodes_count-1)
    {
      /* set slot of last node in snodes to the slot we will place it */
      snodes[snodes_count-1]->bayes_slot = node->parent->bayes_slot;

      /* move this last node to its new slot */
      snodes[node->parent->bayes_slot] = snodes[snodes_count-1];
    }

    /* reset slot of the removed node and decrease count */
    node->parent->bayes_slot = -1;
    --snodes_count;
  }

  /* add select node to the list of speciation nodes */
  node->bayes_slot = snodes_count;
  snodes[snodes_count++] = node;
  node->event = EVENT_SPECIATION;

  /* add left child to coalescent roots unless it is a leaf OR the
     tree rooted at node->left has all branch lengths smaller than minbr */
  if (node->left->edge_count)
  {
    crnodes[crnodes_count] = node->left;
    node->left->bayes_slot = crnodes_count++;
  }

  /* add right child to coalescent roots unless it is a leaf OR the
     tree rooted at node->right has all branch lengths smaller than minbr */
  if (node->right->edge_count)
  {
    crnodes[crnodes_count] = node->right;
    node->right->bayes_slot = crnodes_count++;
  }
}

static void coalesce(unsigned int r)
{
  /*            S                          CR
                *                          *
               / \            ->          / \
              /   \                      /   \
         CR  *     *  CR             C  *     *  C             */

  rtree_t * node = snodes[r];

  /* move the last node of the list to the position of the node
     we just used */
  if (r != (unsigned int)(snodes_count-1))
  {
    snodes[r] = snodes[snodes_count-1];
    snodes[r]->bayes_slot = r;
  }
  --snodes_count;

  /* add the current node to the list of coalescent roots */
  node->bayes_slot = crnodes_count;
  crnodes[crnodes_count++] = node;
  node->event = EVENT_COALESCENT;

  /* remove left child from coalescent roots unless it is a leaf OR the
     tree rooted at node->left has all branch lengths smaller than minbr */
  if (node->left->edge_count)
  {
    /* perform the following only if it is not the last node
       in the list */
    if (node->left->bayes_slot != crnodes_count-1)
    {
      /* set slot of last node in crnodes to the slot we will place it */
      crnodes[crnodes_count-1]->bayes_slot = node->left->bayes_slot;
      
      /* move this last node to its new slot */
      crnodes[node->left->bayes_slot] = crnodes[crnodes_count-1];
    }

    /* reset slot of the removed node and decrease count */
    node->left->bayes_slot = -1;
    crnodes_count--;
  }

  /* now do the same for the right child */
  if (node->right->edge_count)
  {
    /* perform the following only if the parent is not the last node
       in the list */
    if (node->right->bayes_slot != crnodes_count-1)
    {
      /* set slot of last node in crnodes to the slot we will place it */
      crnodes[crnodes_count-1]->bayes_slot = node->right->bayes_slot;
      
      /* move this last node to its new slot */
      crnodes[node->right->bayes_slot] = crnodes[crnodes_count-1];
    }

    /* reset slot of removed node and decrease count */
    node->right->bayes_slot = -1;
    crnodes_count--;
  }

  /* if the parent of the node has two coalescent roots as children
     now, then add it to snodes, i.e. the following case:

              S                               S
              *                               *
             / \                             / \
            /   \                           /   \
       CR  *     *  S          ->      CR  *     *  CR
                / \                             / \
               /   \                           /   \
          CR  *     *  CR                  C  *     *  C      
  */
  if (node->parent && 
      node->parent->left->event == EVENT_COALESCENT && 
      node->parent->right->event == EVENT_COALESCENT)
  {
    assert(node->parent->bayes_slot == -1);

    /* set slot of parent */
    node->parent->bayes_slot = snodes_count;

    /* place parent to the last slot in snodes and increase count */
    snodes[snodes_count++] = node->parent;
  }
}

void bayes(rtree_t * tree,
           int method,
           prior_t * prior,
           struct drand48_data * rstate,
           long seed,
           double * bayes_min_logl,
           double * bayes_max_logl)
{
  int best_index = 0;
  long i;
  long rand_long = 0;
  double rand_double = 0;
  double max = 0;
  double logl = 0;

  *bayes_max_logl = 0;
  *bayes_min_logl = 0;

  if (!opt_quiet)
    fprintf(stdout,"Computing initial delimitation...\n");

  /* check whether all edges are smaller or equal than minbr */
  if (!tree->edge_count)
  {
    fprintf(stderr,"WARNING: All branch lengths are smaller or equal to the "
                   "threshold specified by --minbr. Delimitation equals to "
                   "the null model\n");
    tree->support = 1;
    tree->event = EVENT_COALESCENT;

    return;
  }

  bayes_init(tree, seed);
  
  /* fill DP table */
  dp_recurse(tree, method, prior);

  /* obtain best entry in the root DP table */
  dp_vector_t * vec = tree->vector;
  if (method == PTP_METHOD_MULTI)
  {
    max = vec[0].score_multi;
    for (i = 1; i < tree->edge_count; i++)
    {
      if (max < vec[i].score_multi && vec[i].filled)
      {
        max = vec[i].score_multi;
        best_index = i;
      }
    }
  }
  else
  {
    max = vec[0].score_single;
    for (i = 1; i < tree->edge_count; i++)
    {
      //printf("vec[%d].score_single: %.6f\n", i, vec[i].score_single);
      if (max < vec[i].score_single && vec[i].filled)
      {
        max = vec[i].score_single;
        best_index = i;
      }
    }
  }
  species_count = vec[best_index].species_count;

  unsigned int coal_edge_count = 0;
  unsigned int spec_edge_count;
  double spec_edgelen_sum;
  double coal_edgelen_sum = 0;
  double coal_score = 0;

  if (opt_bayes_startnull && opt_bayes_startrandom)
  {
    fatal("Cannot specify --bayes_startnull and --bayes_startrandom together");
  }
  else if (opt_bayes_startnull)
  {
    tree->event = EVENT_COALESCENT;

    crnodes[crnodes_count++] = tree;
    logl = tree->coal_logl;
    best_index = 0;
    species_count = 1;

    /* set parameters */
    coal_edge_count = tree->edge_count;
    spec_edge_count = 0;
    spec_edgelen_sum = 0;
    coal_edgelen_sum = tree->edgelen_sum;
    coal_score = tree->coal_logl;
    
    /* set all nodes to coalescent */
    init_null(tree);

    /* log log-likelihood at step 0 */
    if (opt_bayes_burnin == 1)
      mcmc_log(logl,species_count);


  }
  else if (opt_bayes_startrandom)
  {
    bool warning_minbr = false;
    logl = random_delimitation(tree,
                               &species_count, 
                               &coal_edge_count, 
                               &coal_edgelen_sum, 
                               &spec_edge_count, 
                               &spec_edgelen_sum, 
                               &coal_score,
                               rstate);
    backtrack_random(tree, &warning_minbr);
    if (warning_minbr)
      fprintf(stderr,"WARNING: A speciation edge is smaller than the specified "
                     "minimum branch length.\n");

    /* log log-likelihood at step 0 */
    if (opt_bayes_burnin == 1)
      mcmc_log(logl,species_count);
  }
  else
  {
    /* ML starting delimitation */
    bool warning_minbr = false;
    backtrack(tree, best_index, &warning_minbr);
    if (warning_minbr)
      fprintf(stderr,"WARNING: A speciation edge is smaller than the specified "
                     "minimum branch length.\n");

    logl = (method == PTP_METHOD_MULTI) ? 
                vec[best_index].score_multi : vec[best_index].score_single;
    
    /* log log-likelihood at step 0 */
    if (opt_bayes_burnin == 1)
      mcmc_log(logl,species_count);
  }

  if (!opt_bayes_startnull && !opt_bayes_startrandom)
  {
    if (method == PTP_METHOD_SINGLE)
    {
      coal_edge_count = tree->edge_count - best_index;
      spec_edge_count = best_index;
      spec_edgelen_sum = tree->vector[best_index].spec_edgelen_sum;
      coal_edgelen_sum = tree->edgelen_sum - spec_edgelen_sum;
    }
    else
    {
      spec_edge_count = best_index;
      spec_edgelen_sum = tree->vector[best_index].spec_edgelen_sum;
      coal_score = tree->vector[best_index].score_multi - 
                        loglikelihood(spec_edge_count, spec_edgelen_sum);
    }
  }

  *bayes_max_logl = logl;
  *bayes_min_logl = logl;

  if (!opt_quiet)
  {
    if (opt_bayes_startnull)
      fprintf(stdout, "Null model log-likelihood: %f\n", logl);
    else if (opt_bayes_startrandom)
      fprintf(stdout, "Random delimitation log-likelihood: %f\n", logl);
    else
      fprintf(stdout, "ML delimitation log-likelihood: %f\n", logl);
  }

  if (opt_bayes_burnin == 1)
  {
    densities[species_count].logl += logl;
  }

  if (opt_bayes_sample == 1)
  {
    if (!opt_quiet)
      printf("1 Log-L: %f\n", logl);
  }

  bayes_stats_init(tree);

  for (i = 1; i < opt_bayes_runs; ++i)
  {

    //if (opt_bayes_burnin == i)
    //  bayes_stats_init(tree);

    /* throw a coin to decide whether to convert a coalescent root to a
       speciation or the other way round */
    drand48_r(rstate, &rand_double);
    int speciation = (rand_double >= 0.5) ? 1 : 0;

    if ((speciation && crnodes_count) || (snodes_count == 0))
    {

      /*            CR                         S
                    *                          *
                   / \            ->          / \
                  /   \                      /   \
              C  *     *  C             CR  *     *  CR            */


      /* select a coalescent root, split it into two coalescent nodes */
      lrand48_r(rstate, &rand_long);
      unsigned int r = (unsigned int)(rand_long % crnodes_count);
      rtree_t * node = crnodes[r];

      /* store the count of crnodes for the Hasting ratio */
      double old_crnodes_count = crnodes_count;

      /* speciate */
      speciate(r);

      /* store the new count of snodes for the Hasting ratio */
      double new_snodes_count = snodes_count;

      /* TODO: distinguish between single- and multi-rate methods */

      /* subtract the two edges (left and right) from the coalescent
         distribution and add them to the speciation distribution */
      int edge_count_diff = 0;
      double edgelen_sum_diff = 0;
      if (node->left->length > opt_minbr)
      {
        ++edge_count_diff;
        edgelen_sum_diff += node->left->length;
      }

      if (node->right->length > opt_minbr)
      {
        ++edge_count_diff;
        edgelen_sum_diff += node->right->length;
      }
      
      if (method == PTP_METHOD_SINGLE)
      {
        coal_edgelen_sum -= edgelen_sum_diff;
        coal_edge_count -= edge_count_diff;
      }
      spec_edgelen_sum += edgelen_sum_diff;
      spec_edge_count += edge_count_diff;

      /* compute new log-likelihood */
      double new_logl;
      if (spec_edge_count == 0 || (method == PTP_METHOD_SINGLE && coal_edge_count == 0))
        new_logl = tree->coal_logl;
      else
      {
        assert((method == PTP_METHOD_MULTI) || (coal_edge_count > 0));
        assert(spec_edge_count > 0);
        if (method == PTP_METHOD_SINGLE)
          new_logl = loglikelihood(coal_edge_count, coal_edgelen_sum) +
                     loglikelihood(spec_edge_count, spec_edgelen_sum);
        else
          new_logl = coal_score - node->coal_logl + 
                     node->left->coal_logl + node->right->coal_logl +
                     loglikelihood(spec_edge_count, spec_edgelen_sum);
          
      }

      if (new_logl > *bayes_max_logl)
        *bayes_max_logl = new_logl;
      if (i+1 < opt_bayes_burnin)
        *bayes_min_logl = *bayes_max_logl;
      else if (new_logl < *bayes_min_logl)
        *bayes_min_logl = new_logl;

      /* Hastings ratio */
      double a = exp(new_logl - logl) * (old_crnodes_count / new_snodes_count);

      /* update frequencies */
      if (i+1 >= opt_bayes_burnin)
      {
        densities[species_count+1].logl += new_logl;
      }

      /* decide whether to accept or reject proposal */
      drand48_r(rstate, &rand_double);
      if (rand_double <= a)
      {
        /* accept */
        if ((i+1) % opt_bayes_sample == 0)
        {
          if (!opt_quiet)
            printf("%ld Log-L: %f\n", i+1, new_logl);
          if (i+1 >= opt_bayes_burnin)
            mcmc_log(new_logl,species_count+1);
        }

        /* update support values information */
        if (i+1 >= opt_bayes_burnin)
          node->speciation_start = i;

        accept_count++;
        species_count++;
        logl = new_logl;
        if (method == PTP_METHOD_MULTI)
          coal_score = coal_score - node->coal_logl +
                       node->left->coal_logl + node->right->coal_logl;
        continue;
      }
      else
      { 
        /* reject */
        if ((i+1) % opt_bayes_sample == 0)
        {
          if (!opt_quiet)
            printf("%ld Log-L: %f\n", i+1, new_logl);
          if (i+1 >= opt_bayes_burnin)
            mcmc_log(new_logl,species_count+1);
        }

        if (i+1 >= opt_bayes_burnin)
          node->speciation_count++;

        if (method == PTP_METHOD_SINGLE)
        {
          coal_edgelen_sum += edgelen_sum_diff;
          coal_edge_count += edge_count_diff;
        }
        spec_edgelen_sum -= edgelen_sum_diff;
        spec_edge_count -= edge_count_diff;
        coalesce(node->bayes_slot);
      }
    }
    else
    {

      /*            S                          CR
                    *                          *
                   / \            ->          / \
                  /   \                      /   \
             CR  *     *  CR             C  *     *  C         */

      lrand48_r(rstate, &rand_long);
      unsigned int r = (unsigned int)(rand_long % snodes_count);
      rtree_t * node = snodes[r];

      /* store the count of snodes for the Hastings ratio */
      double old_snodes_count = snodes_count;

      /* coalesce */
      coalesce(r);

      double new_crnodes_count = crnodes_count;

      /* TODO: distinguish between single- and multi-rate methods */

      /* subtract the two edges (left and right) from the speciation
         distribution and add them to the coalescent distribution */
      int edge_count_diff = 0;
      double edgelen_sum_diff = 0;
      if (node->left->length > opt_minbr)
      {
        ++edge_count_diff;
        edgelen_sum_diff += node->left->length;
      }

      if (node->right->length > opt_minbr)
      {
        ++edge_count_diff;
        edgelen_sum_diff += node->right->length;
      }
      if (method == PTP_METHOD_SINGLE)
      {
        coal_edgelen_sum += edgelen_sum_diff;
        coal_edge_count += edge_count_diff;
      }
      spec_edgelen_sum -= edgelen_sum_diff;
      spec_edge_count -= edge_count_diff;

      /* compute new log-likelihood */
      double new_logl;
      if (spec_edge_count == 0 || (method == PTP_METHOD_SINGLE && coal_edge_count == 0))
        new_logl = tree->coal_logl;
      else
      {
        assert((method == PTP_METHOD_MULTI) || (coal_edge_count > 0));
        assert(spec_edge_count > 0);
        if (method == PTP_METHOD_SINGLE)
          new_logl = loglikelihood(coal_edge_count, coal_edgelen_sum) +
                     loglikelihood(spec_edge_count, spec_edgelen_sum);
        else
          new_logl = coal_score - node->left->coal_logl - node->right->coal_logl + 
                     node->coal_logl +
                     loglikelihood(spec_edge_count, spec_edgelen_sum);

      }

      if (new_logl > *bayes_max_logl)
        *bayes_max_logl = new_logl;
      if (i+1 < opt_bayes_burnin)
        *bayes_min_logl = *bayes_max_logl;
      else if (new_logl < *bayes_min_logl)
        *bayes_min_logl = new_logl;

      /* Hastings ratio */
      double a = exp(new_logl - logl) * (old_snodes_count / new_crnodes_count);

      /* update frequencies */
      if (i+1 >= opt_bayes_burnin)
      {
        densities[species_count-1].logl += new_logl;
      }

      /* decide whether to accept or reject proposal */
      drand48_r(rstate, &rand_double);
      if (rand_double <= a)
      {
        /* accept */
        if ((i+1) % opt_bayes_sample == 0)
        {
          if (!opt_quiet)
            printf("%ld Log-L: %f\n", i+1, new_logl);
          if (i+1 >= opt_bayes_burnin)
            mcmc_log(new_logl,species_count-1);
        }

        /* update support values information */
        if (i+1 >= opt_bayes_burnin)
        {
          node->speciation_count = node->speciation_count + 
                                   i - node->speciation_start;
          node->speciation_start = -1;
        }

        accept_count++;
        species_count--;
        logl = new_logl;
        if (method == PTP_METHOD_MULTI)
          coal_score = coal_score - node->left->coal_logl - node->right->coal_logl + 
                       node->coal_logl;

        continue;
      }
      else
      {
        /* reject */
        if ((i+1) % opt_bayes_sample == 0)
        {
          if (!opt_quiet)
            printf("%ld Log-L: %f\n", i+1, new_logl);
          if (i+1 >= opt_bayes_burnin)
            mcmc_log(new_logl,species_count-1);
        }
        if (method == PTP_METHOD_SINGLE)
        {
          coal_edgelen_sum -= edgelen_sum_diff;
          coal_edge_count -= edge_count_diff;
        }
        spec_edgelen_sum += edgelen_sum_diff;
        spec_edge_count += edge_count_diff;
        speciate(node->bayes_slot);
        if (i+1 >= opt_bayes_burnin)
        {
          node->speciation_count--;
        }
      }
    }
  }

  /* TODO: DEBUG variables for checking the max likelihood bayesian runs give.
     Must be removed */
  bayes_finalize(tree, *bayes_min_logl, *bayes_max_logl, seed);

}
