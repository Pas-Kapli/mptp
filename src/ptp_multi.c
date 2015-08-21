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
#include <gsl/gsl_cdf.h>

double loglikelihood(int edge_count, double edgelen_sum)
{
  assert(edge_count >= 0);

  if (edge_count == 0 || edgelen_sum == 0) // TODO: Try to find out if sum == 0 makes sense here
    return 0;

  return edge_count * (log(edge_count) - 1 - log(edgelen_sum));
}

void init_tree_data(rtree_t * tree)
{

  /* initialize everything to coalescent, speciations will then
     over-write this value */
  tree->event = EVENT_COALESCENT;

  if (tree->left)
    init_tree_data(tree->left);
  if (tree->right)
    init_tree_data(tree->right);

  node_information* info = malloc(sizeof(node_information));

  // TODO: Check whether this is the best way to handle those
  //   nasty zero-length edges.

  info->spec_array = calloc(tree->valid_edge_count + 1, sizeof(spec_entry));

  int i;
  for (i = 0; i <= tree->valid_edge_count; i++)
  {
    info->spec_array[i].taken_left_index = -1;
    info->spec_array[i].taken_right_index = -1;
    info->spec_array[i].score_multi = -999999999;
    info->spec_array[i].score_single = -999999999;
    info->spec_array[i].sum_speciation_edges_subtree = 0;
    info->spec_array[i].coalescent_value = 0;
    info->spec_array[i].speciation_value = 0;
    info->spec_array[i].num_species = 0;
    info->spec_array[i].valid = false;
  }

  assert(tree->valid_edge_count >= 0);

  tree->coalescent_logl = loglikelihood(tree->valid_edge_count,
                                        tree->valid_edgelen_sum);
  tree->data = info;
}

void free_tree_data(rtree_t * tree)
{
  if (tree->left)
  {
    free_tree_data(tree->left);
  }
  if (tree->right)
  {
    free_tree_data(tree->right);
  }
  free(((node_information*) (tree->data))->spec_array);
  free((node_information*) (tree->data));
  tree->data = NULL;
}

static void init_additional_tree_data(rtree_t * tree)
{
  // pre-order traversal
  int num_known_speciation_edges = 0;
  double sum_known_speciation_edges = 0;
  if (tree->parent)
  {
    node_information* parent_data = (node_information*) (tree->parent->data);

    if (tree->length <= opt_minbr)
    {
      num_known_speciation_edges = parent_data->num_known_speciation_edges;
      sum_known_speciation_edges = parent_data->sum_known_speciation_edges;
    }
    else
    {
      num_known_speciation_edges = parent_data->num_known_speciation_edges + 1;
      sum_known_speciation_edges =
        parent_data->sum_known_speciation_edges + tree->length;
    }

    if (tree == tree->parent->left)
    {
        if (tree->parent->right && tree->parent->right->length > opt_minbr)
        {
          num_known_speciation_edges++;
          sum_known_speciation_edges += tree->parent->right->length;
        }
    }
    else
    {
      if (tree->parent->left && tree->parent->left->length > opt_minbr)
      {
        num_known_speciation_edges++;
        sum_known_speciation_edges += tree->parent->left->length;
      }
    }
  }

  ((node_information*) (tree->data))->num_known_speciation_edges =
    num_known_speciation_edges;
  ((node_information*) (tree->data))->sum_known_speciation_edges =
    sum_known_speciation_edges;

  if (tree->left)
  {
    init_additional_tree_data(tree->left);
  }
  if (tree->right)
  {
    init_additional_tree_data(tree->right);
  }
}

void multi_traversal(rtree_t * tree, bool multiple_lambda,
  PRIOR_FUNC species_logprior, prior_inf prior_information)
{
  // post order traversal
  if (tree->left)
  {
    multi_traversal(tree->left, multiple_lambda, species_logprior, prior_information);
  }
  if (tree->right)
  {
    multi_traversal(tree->right, multiple_lambda, species_logprior, prior_information);
  }

  node_information* data = (node_information*) (tree->data);
  spec_entry* spec_array_act = data->spec_array;

  // zero speciation edges means that the current node is the MRCA of a species
  spec_array_act[0].coalescent_value = tree->coalescent_logl;
  double combined_spec_sum = data->sum_known_speciation_edges;
  int combined_spec_num = data->num_known_speciation_edges;
  assert(combined_spec_num >= 0);
  double speciation_value =
    loglikelihood(combined_spec_num, combined_spec_sum)
    + species_logprior(1, prior_information);
  spec_array_act[0].speciation_value = speciation_value;
  spec_array_act[0].sum_speciation_edges_subtree = 0;
  spec_array_act[0].score_multi = tree->coalescent_logl + speciation_value;
  spec_array_act[0].score_single = tree->coalescent_logl + speciation_value;
  spec_array_act[0].num_species = 1;
  spec_array_act[0].valid = true;

  // now we fill the arrays for the current node being a speciation event
  if (tree->left && tree->right)
  // inner node with two children
  {
    spec_entry* spec_array_left =
      ((node_information*) (tree->left->data))->spec_array;
    spec_entry* spec_array_right =
      ((node_information*) (tree->right->data))->spec_array;

    int num_valid_child_edges = 2;
    if (tree->left->length <= opt_minbr)
    {
      num_valid_child_edges--;
    }
    if (tree->right->length <= opt_minbr)
    {
      num_valid_child_edges--;
    }

    int i;
    for (i = 0; i <= tree->left->valid_edge_count; i++)
    {
      int j;
      for (j = 0; j <= tree->right->valid_edge_count; j++)
      {
        if (spec_array_left[i].valid && spec_array_right[j].valid)
        {
          int current_num_species = spec_array_left[i].num_species
            + spec_array_right[j].num_species;

          double combined_spec_sum = 0;
          int combined_spec_num = 0;
          double sum_speciation_edges_subtree = 0;

          if (tree->left->length > opt_minbr)
          {
            combined_spec_sum += tree->left->length;
            combined_spec_num++;
            sum_speciation_edges_subtree += tree->left->length;
          }
          if (tree->right->length > opt_minbr)
          {
            combined_spec_sum += tree->right->length;
            combined_spec_num++;
            sum_speciation_edges_subtree += tree->right->length;
          }

          combined_spec_sum += data->sum_known_speciation_edges
            + spec_array_left[i].sum_speciation_edges_subtree
            + spec_array_right[j].sum_speciation_edges_subtree;
          combined_spec_num += data->num_known_speciation_edges
            + i + j;

          sum_speciation_edges_subtree +=
            spec_array_left[i].sum_speciation_edges_subtree
            + spec_array_right[j].sum_speciation_edges_subtree;

          double sum_coalescent_edges = tree->valid_edgelen_sum
            - sum_speciation_edges_subtree;
          double num_coalescent_edges =
            tree->valid_edge_count - (i + j + num_valid_child_edges);

          double coalescent_value_multi = spec_array_left[i].coalescent_value
              + spec_array_right[j].coalescent_value;

          if (num_coalescent_edges < 0)
            printf("Error: %d %d %d %d\n", tree->valid_edge_count, i, j, num_valid_child_edges);
          assert(num_coalescent_edges >= 0);
          double coalescent_value_single =
            loglikelihood(num_coalescent_edges, sum_coalescent_edges);

          assert(combined_spec_num >= 0);
          double speciation_value =
            loglikelihood(combined_spec_num, combined_spec_sum)
            + species_logprior(current_num_species, prior_information);

          double score_multi = coalescent_value_multi + speciation_value;
          double score_single = coalescent_value_single + speciation_value;

          double score = score_multi;
          double old_score = spec_array_act[i+j+num_valid_child_edges].score_multi;
          if (!multiple_lambda)
          {
            score = score_single;
            old_score = spec_array_act[i+j+num_valid_child_edges].score_single;
          }

          if (score > old_score)
          {
            spec_array_act[i+j+num_valid_child_edges].score_multi = score_multi;
            spec_array_act[i+j+num_valid_child_edges].score_single = score_single;
            spec_array_act[i+j+num_valid_child_edges].sum_speciation_edges_subtree =
              sum_speciation_edges_subtree;
            spec_array_act[i+j+num_valid_child_edges].coalescent_value =
              coalescent_value_multi;
            spec_array_act[i+j+num_valid_child_edges].speciation_value =
              speciation_value;
            spec_array_act[i+j+num_valid_child_edges].taken_left_index = i;
            spec_array_act[i+j+num_valid_child_edges].taken_right_index = j;
            spec_array_act[i+j+num_valid_child_edges].num_species =
              current_num_species;
            spec_array_act[i+j+num_valid_child_edges].valid = true;
          }
        }
      }
    }
  }
}

void print_subtree_leaves(rtree_t * tree)
{
  if (tree->left)
  {
    print_subtree_leaves(tree->left);
  }
  if (tree->right)
  {
    print_subtree_leaves(tree->right);
  }

  if (!tree->left && !tree->right)
  {
    printf("%s\n",tree->label);
  }
}

int current_species_idx = 1;

static void backtrack_species_assignment(rtree_t * tree, int pos,
  bool *print_speciation_warning)
{
  spec_entry* spec_array = ((node_information*) (tree->data))->spec_array;
  if (spec_array[pos].taken_left_index != -1
    && spec_array[pos].taken_right_index != -1)
    // speciation event
  {
    tree->event = EVENT_SPECIATION;
    if (tree->length <= opt_minbr && tree->parent)
    {
      (*print_speciation_warning) = true;
    }
    backtrack_species_assignment(tree->left, spec_array[pos].taken_left_index,
      print_speciation_warning);
    backtrack_species_assignment(tree->right,
      spec_array[pos].taken_right_index, print_speciation_warning);
  }
  else  // coalescent event
  {
    tree->event = EVENT_COALESCENT;
    if (opt_quiet)
    {
      current_species_idx++;
    }
    else
    {
      printf("\nSpecies %d:\n", current_species_idx++);
      print_subtree_leaves(tree);
    }
  }
}

void print_null_model(rtree_t * tree)
{
  printf("Species 1:\n");
  current_species_idx++;
  print_subtree_leaves(tree);
}

bool likelihood_ratio_test(rtree_t * tree, double likelihood_alternative_model,
  double p_value, int degrees_of_freedom, double * computed_p_value)
{
  double likelihood_null_model = tree->coalescent_logl;
  double lr = 2 * (likelihood_alternative_model - likelihood_null_model);
  (*computed_p_value) = 1 - gsl_cdf_chisq_P(lr, degrees_of_freedom);
    // should return
    //  1/(2**(df/2) * gamma(df/2)) * integral(t**(df/2-1) * exp(-t/2), t=0..x)
    // see http://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chdtr.html#scipy.special.chdtr

  bool good_delimitation = true;
  if ((*computed_p_value) > p_value)
  {
    good_delimitation = false;
  }
  return good_delimitation;
}

delimit_stats * ptp_multi_heuristic(rtree_t * tree,
                                    bool multiple_lambda,
                                    PRIOR_FUNC prior_function,
                                    prior_inf prior_information)
{
  init_tree_data(tree);
  init_additional_tree_data(tree);
  multi_traversal(tree, multiple_lambda, prior_function, prior_information);
  //print_debug_information(tree);

  double max = 0;
  node_information* data = (node_information*) (tree->data);
  spec_entry* spec_array = data->spec_array;
  int i;
  int pos = 0;

  if (multiple_lambda)
  {
    for (i = 0; i < tree->valid_edge_count; i++)
    {
      if (max < spec_array[i].score_multi)
      {
        max = spec_array[i].score_multi;
        pos = i;
      }
    }
  }
  else // single lambda
  {
    for (i = 0; i < tree->valid_edge_count; i++)
    {
      if (max < spec_array[i].score_single)
      {
        max = spec_array[i].score_single;
        pos = i;
      }
    }
  }

  if (!opt_quiet)
  {
    printf("Num big enough edges in tree: %d\n",
      tree->valid_edge_count);
    printf("Score Null Model: %.6f\n", tree->coalescent_logl);
    printf("Best score found single: %.6f\n", spec_array[pos].score_single);
    printf("Best score found multi: %.6f\n", spec_array[pos].score_multi);
  }

  bool good_delimitation = true;
  double computed_p_value = -1;
  if (multiple_lambda)
  {
    good_delimitation = likelihood_ratio_test(tree, spec_array[pos].score_multi,
      opt_pvalue, 1, &computed_p_value);
    // TODO: Find out whether here is really just 1 degree of freedom or not
  }
  else
  {
    good_delimitation = likelihood_ratio_test(tree,
      spec_array[pos].score_single, opt_pvalue, 1, &computed_p_value);
  }

  if (!opt_quiet)
  {
    printf("Computed P-value: %.6f\n", computed_p_value);
  }

  if (good_delimitation)
  {
    bool print_speciation_warning = false;
    backtrack_species_assignment(tree, pos,
      &print_speciation_warning);
    if (print_speciation_warning && !opt_quiet)
    {
      printf("WARNING: A speciation edge in the result is too small.\n");
    }
  }
  else
  {
    if (!opt_quiet)
    {
      printf("The Null Model is the preferred one.\n");
      print_null_model(tree);
    }
    else
    {
      current_species_idx = 2;
    }
  }
  if (!opt_quiet)
  {
    printf("Number of delimited species: %d\n", current_species_idx - 1);
  }
  if (tree->valid_edge_count == 0 && !opt_quiet)
  {
    fprintf(stderr, "WARNING: The tree has no edges > %f. All edges have been ignored. \n", opt_minbr);
  }

  delimit_stats* solution = malloc(sizeof(delimit_stats));
  solution->number_of_species = current_species_idx - 1;
  solution->is_null_model = !good_delimitation;
  if (good_delimitation)
  {
    solution->score_single = spec_array[pos].score_single;
    solution->score_multi = spec_array[pos].score_multi;
  }
  else
  {
    solution->score_single = tree->coalescent_logl;
    solution->score_multi = tree->coalescent_logl;
  }

  free_tree_data(tree);
  return solution;
}
