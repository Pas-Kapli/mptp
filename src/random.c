#include "delimit.h"

static long min_species;
static long max_species;
static long species_count;

static int cb_node_select(rtree_t * node)
{
  if (!node->edge_count) return 0;

  /* check if not selecting node is possible */
  if (min_species+1 > species_count)
  {
    /* we must select the node */
    node->event = EVENT_COALESCENT;
    max_species = max_species - node->max_species_count + 1;
    return 0;
  }
  
  /* check if selecting the node is possible */
  if (max_species - node->max_species_count + 1 < species_count)
  {
    /* we must NOT select the node */
    node->event = EVENT_SPECIATION;
    min_species = min_species+1;
    return 1;
  }

  /* otherwise, we just throw a coin and select one of the two cases */
  if (drand48() >= 0.5)
  {
    /* don't select */
    node->event = EVENT_SPECIATION;
    min_species = min_species+1;
    return 1;
  }

  /* otherwise select node */
  node->event = EVENT_COALESCENT;
  max_species = max_species - node->max_species_count + 1;
  return 0;
}

double random_delimitation(rtree_t * root,
                           long * delimited_species,
                           unsigned int * coal_edge_count,
                           double * coal_edgelen_sum,
                           unsigned int * spec_edge_count,
                           double * spec_edgelen_sum,
                           double * coal_score)
{
  long i;
  double logl = 0;
  int edge_count = 0;
  double edgelen_sum = 0;

  /* initialize */
  min_species = 1;
  max_species = root->max_species_count;


  species_count = (lrand48() % root->max_species_count) + 1;

  rtree_t ** inner_node_list =  (rtree_t **)xmalloc(species_count *
                                                    sizeof(rtree_t *));

  long count = rtree_traverse(root, cb_node_select, inner_node_list);

  for (i = 0; i < count; ++i)
  {
    logl += inner_node_list[i]->coal_logl;
    edge_count += inner_node_list[i]->edge_count;
    edgelen_sum += inner_node_list[i]->edgelen_sum;
  }
  *coal_score = logl;

  /* if we have PTP single logl is different */
  if (opt_ml_single || opt_bayes_single)
    logl = loglikelihood(edge_count, edgelen_sum);
  
  
  /* append speciation part log-likelihood */
  logl += loglikelihood(root->edge_count - edge_count,
                        root->edgelen_sum - edgelen_sum);

  free(inner_node_list);

  assert(count == species_count);

  *delimited_species = species_count;
  *coal_edge_count = edge_count;
  *coal_edgelen_sum = edgelen_sum;
  *spec_edge_count = root->edge_count - edge_count;
  *spec_edgelen_sum = root->edgelen_sum - edgelen_sum;

  return logl;
}


