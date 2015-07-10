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

double compute_loglikelihood(int num, double sum)
{
  return num * log(num) - num - num * log(sum);
}

void init_tree_data(rtree_t * tree)
{
  // post-order traversal
  int subtree_size_edges = 0;
  double subtree_sum_edges = 0;
  if (tree->left)
  {
    init_tree_data(tree->left);
    node_information* left_data = (node_information*) (tree->left->data);
    subtree_size_edges += left_data->num_edges_subtree + 1;
    subtree_sum_edges += left_data->sum_edges_subtree;
    subtree_sum_edges += tree->left->length;
  }
  if (tree->right)
  {
    init_tree_data(tree->right);
    node_information* right_data = (node_information*) (tree->right->data);
    subtree_size_edges += right_data->num_edges_subtree + 1;
    subtree_sum_edges += right_data->sum_edges_subtree;
    subtree_sum_edges += tree->right->length;
  }
  node_information* info = malloc(sizeof(node_information));
  info->num_edges_subtree = subtree_size_edges;
  info->sum_edges_subtree = subtree_sum_edges;
  info->spec_array = calloc(subtree_size_edges + 1, sizeof(spec_entry));

  int i;
  for (i = 0; i <= subtree_size_edges; i++)
  {
    info->spec_array[i].taken_left_index = -1;
    info->spec_array[i].taken_right_index = -1;
    info->spec_array[i].score_multi = -999999999;
    info->spec_array[i].score_single = -999999999;
    info->spec_array[i].sum_speciation_edges_subtree = 0;
    info->spec_array[i].coalescent_value = 0;
    info->spec_array[i].speciation_value = 0;
  }

  info->coalescent =
    compute_loglikelihood(subtree_size_edges,subtree_sum_edges);
  tree->data = info;
}

void init_additional_tree_data(rtree_t * tree)
{
  // pre-order traversal
  int num_known_speciation_edges = 0;
  double sum_known_speciation_edges = 0;
  if (tree->parent)
  {
    num_known_speciation_edges =
      ((node_information*) (tree->parent->data))->num_known_speciation_edges
        + 1;
    sum_known_speciation_edges =
      ((node_information*) (tree->parent->data))->sum_known_speciation_edges
        + tree->length;

    if (tree == tree->parent->left)
    {
      if (tree->parent->right)
      {
        num_known_speciation_edges++;
        sum_known_speciation_edges += tree->parent->right->length;
      }
    }
    else
    {
      if (tree->parent->left)
      {
        num_known_speciation_edges++;
        sum_known_speciation_edges += tree->parent->left->length;
      }
    }
  }
  if (tree->left)
  {
    num_known_speciation_edges++;
    sum_known_speciation_edges += tree->left->length;
  }
  if (tree->right)
  {
    num_known_speciation_edges++;
    sum_known_speciation_edges += tree->right->length;
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

void multi_traversal(rtree_t * tree, bool multiple_lambda)
{
  // post order traversal
  if (tree->left)
  {
    multi_traversal(tree->left, multiple_lambda);
  }
  if (tree->right)
  {
    multi_traversal(tree->right, multiple_lambda);
  }

  if (tree->left && tree->right)
  {
    node_information* data = (node_information*) (tree->data);
    node_information* data_left = (node_information*) (tree->left->data);
    node_information* data_right = (node_information*) (tree->right->data);
    spec_entry* spec_array_act = data->spec_array;
    spec_entry* spec_array_left =
      ((node_information*) (tree->left->data))->spec_array;
    spec_entry* spec_array_right =
      ((node_information*) (tree->right->data))->spec_array;

    spec_array_act[0].coalescent_value = data->coalescent;
    spec_array_act[0].score_multi = data->coalescent;
    spec_array_act[0].score_single = data->coalescent;
    spec_array_act[0].speciation_value = 0;
    spec_array_act[0].sum_speciation_edges_subtree = 0;

    int i;
    for (i = 0; i <= data_left->num_edges_subtree; i+=2)
    {
      if (spec_array_left[i].score_multi == -999999999)
      {
        break;
      }
      int j;
      for (j = 0; j <= data_right->num_edges_subtree; j+=2)
      {
        if (spec_array_right[j].score_multi == -999999999)
        {
          break;
        }
        double combined_spec_sum = data->sum_known_speciation_edges
          + spec_array_left[i].sum_speciation_edges_subtree
          + spec_array_right[j].sum_speciation_edges_subtree;
        int combined_spec_num = data->num_known_speciation_edges
          + i + j;

        double sum_speciation_edges_subtree =
          spec_array_left[i].sum_speciation_edges_subtree
          + spec_array_right[j].sum_speciation_edges_subtree
          + tree->left->length
          + tree->right->length;

        double sum_coalescent_edges = data->sum_edges_subtree
          - sum_speciation_edges_subtree;
        double num_coalescent_edges =
          data->num_edges_subtree - (i+j+2);

        double coalescent_value_multi = spec_array_left[i].coalescent_value
            + spec_array_right[j].coalescent_value;

        double coalescent_value_single = compute_loglikelihood(num_coalescent_edges,
            sum_coalescent_edges);

        double speciation_value = compute_loglikelihood(combined_spec_num, combined_spec_sum);

        double score_multi = coalescent_value_multi + speciation_value;
        double score_single = coalescent_value_single + speciation_value;

        double score = score_multi;
        double old_score = spec_array_act[i+j+2].score_multi;
        if (!multiple_lambda)
        {
          score = score_single;
          old_score = spec_array_act[i+j+2].score_single;
        }

        if (score > old_score)
        {
          spec_array_act[i+j+2].score_multi = score_multi;
          spec_array_act[i+j+2].score_single = score_single;
          spec_array_act[i+j+2].sum_speciation_edges_subtree =
            sum_speciation_edges_subtree;
          spec_array_act[i+j+2].coalescent_value = coalescent_value_multi;
          spec_array_act[i+j+2].speciation_value = speciation_value;
          spec_array_act[i+j+2].taken_left_index = i;
          spec_array_act[i+j+2].taken_right_index = j;
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

int current_species_num = 1;

void backtrack_species_assignment(rtree_t * tree, int pos)
{
  spec_entry* spec_array = ((node_information*) (tree->data))->spec_array;
  if (spec_array[pos].taken_left_index != -1
    && spec_array[pos].taken_right_index != -1)
    // speciation event
  {
    backtrack_species_assignment(tree->left, spec_array[pos].taken_left_index);
    backtrack_species_assignment(tree->right,
      spec_array[pos].taken_right_index);
  }
  else  // coalescent event
  {
    printf("\nSpecies %d:\n", current_species_num++);
    print_subtree_leaves(tree);
  }
}

void ptp_multi_heuristic(rtree_t * tree, bool multiple_lambda)
{
  init_tree_data(tree);
  init_additional_tree_data(tree);
  multi_traversal(tree, multiple_lambda);
  double max = 0;
  node_information* data = (node_information*) (tree->data);
  spec_entry* spec_array = data->spec_array;
  int i;
  int pos = 0;
  //printf("Debug: Whole score array\n");
  for (i = 0; i < data->num_edges_subtree; i++)
  {
    //printf("%.6f ", spec_array[i].score);
    if (max < spec_array[i].score_single)
    {
      max = spec_array[i].score_single;
      pos = i;
    }
  }
  //printf("\n");

  printf("Null logl: %.6f\n", data->coalescent);
  printf("Best score found single: %.6f\n", spec_array[pos].score_single);
  printf("Best score found multi: %.6f\n", spec_array[pos].score_multi);

  /*printf("Debug position: %d\n", pos);

  printf("Debug sum_speciation_edges_subtree: %.6f\n", spec_array[pos].sum_speciation_edges_subtree);
  printf("Debug sum all edges: %.6f\n", data->sum_edges_subtree);
  printf("Debug sum_coalescent_edges_subtree: %.6f\n", data->sum_edges_subtree - spec_array[pos].sum_speciation_edges_subtree);
  printf("Debug coalescent value: %.6f\n", spec_array[pos].coalescent_value);
  printf("Debug speciation value: %.6f\n", spec_array[pos].speciation_value);
  printf("Debug sum known speciation edges: %.6f\n", data->sum_known_speciation_edges);
  printf("Debug num known speciation edges: %d\n", data->num_known_speciation_edges);

  spec_entry* spec_array_left = ((node_information*) (tree->left->data))->spec_array;
  spec_entry* spec_array_right = ((node_information*) (tree->right->data))->spec_array;
  i = spec_array[pos].taken_left_index;
  int j = spec_array[pos].taken_right_index;

  double combined_spec_sum = data->sum_known_speciation_edges
    + spec_array_left[i].sum_speciation_edges_subtree
    + spec_array_right[j].sum_speciation_edges_subtree;
  int combined_spec_num = data->num_known_speciation_edges + i + j;

  printf("Debug taken left index: %d\n", spec_array[pos].taken_left_index);
  printf("Debug taken right index: %d\n", spec_array[pos].taken_right_index);

  printf("Debug combined spec sum: %.6f\n", combined_spec_sum);
  printf("Debug combined spec num: %d\n", combined_spec_num);

  printf("Strange score right: %.6f\n", spec_array_right[j].score);
  printf("Strange sum_speciation_edges_subtree right: %.6f\n", spec_array_right[j].sum_speciation_edges_subtree);
  printf("Strange speciation value right: %.6f\n", spec_array_right[j].speciation_value);
  printf("Strange coalescent value right: %.6f\n", spec_array_right[j].coalescent_value);

  printf("Debug: Whole score array RIGHT\n");
  for (i = 0; i < ((node_information*) (tree->right->data))->num_edges_subtree; i++)
  {
    printf("%.6f ", spec_array_right[i].score);
  }
  printf("\n");*/


  backtrack_species_assignment(tree, pos);
}
