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

static long min_species;
static long max_species;
static long species_count;
static unsigned short * g_rstate;

static int cb_node_select(rtree_t * node)
{
  double rand_double = 0;

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
  rand_double = mptp_erand48(g_rstate);
  if (rand_double >= 0.5)
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
                           long * coal_edge_count,
                           double * coal_edgelen_sum,
                           long * spec_edge_count,
                           double * spec_edgelen_sum,
                           double * coal_score,
                           unsigned short * rstate)
{
  int edge_count = 0;
  long i;
  long rand_long = 0;
  double logl = 0;
  double edgelen_sum = 0;

  /* initialize */
  min_species = 1;
  max_species = root->max_species_count;
  g_rstate = rstate;

  rand_long = mptp_nrand48(rstate);
  if (!root->max_species_count)
    species_count = (rand_long % root->leaves) + 1;
  else
    species_count = (rand_long % root->max_species_count) + 1;

  rtree_t ** inner_node_list =  (rtree_t **)xmalloc((size_t)species_count *
                                                    sizeof(rtree_t *));

  long count = rtree_traverse(root, cb_node_select, rstate, inner_node_list);

  for (i = 0; i < count; ++i)
  {
    logl += inner_node_list[i]->coal_logl;
    edge_count += inner_node_list[i]->edge_count;
    edgelen_sum += inner_node_list[i]->edgelen_sum;
  }
  *coal_score = logl;

  /* if we have PTP single logl is different */
  if (opt_method == PTP_METHOD_SINGLE)
    logl = loglikelihood(edge_count, edgelen_sum);


  /* append speciation part log-likelihood */
  logl += loglikelihood(root->edge_count - edge_count,
                        root->edgelen_sum - edgelen_sum);

  free(inner_node_list);
  
  assert(count <= species_count);
  if (count < species_count)
  {
    /* TODO: This fixes issue #82, but we should implement a better, non-biased
       way of generatng random starting delimitations */
    species_count = count;
  }
  assert(count == species_count);

  *delimited_species = species_count;
  *coal_edge_count = edge_count;
  *coal_edgelen_sum = edgelen_sum;
  *spec_edge_count = root->edge_count - edge_count;
  *spec_edgelen_sum = root->edgelen_sum - edgelen_sum;

  return logl;
}
