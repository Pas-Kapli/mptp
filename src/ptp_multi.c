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

  info->coalescent = subtree_size_edges * log(subtree_size_edges)
                      - subtree_size_edges
                      - subtree_size_edges * log(subtree_sum_edges);
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
      ((node_information*) (tree->parent->data))->num_known_speciation_edges + 1;
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

double compute_speciation_value(int num, double sum)
{
  return num * log(num) - num - num * log(sum);
}

void multi_traversal(rtree_t * tree)
{
  // post order traversal
  if (tree->left)
  {
    multi_traversal(tree->left);
  }
  if (tree->right)
  {
    multi_traversal(tree->right);
  }

  if (tree->left && tree->right)
  {
    node_information* data = (node_information*) (tree->data);
    spec_entry* spec_array_act = data->spec_array;
    spec_entry* spec_array_left =
      ((node_information*) (tree->left->data))->spec_array;
    spec_entry* spec_array_right =
      ((node_information*) (tree->right->data))->spec_array;

    spec_array_act[0].coalescent_sum_subtree = data->coalescent;

    int i;
    for (i = 0; i <= ((node_information*) tree->left->data)->num_edges_subtree; i++)
    {
      int j;
      for (j = 0; j <= ((node_information*) tree->right->data)->num_edges_subtree; j++)
      {
        double combined_spec_sum = data->sum_known_speciation_edges
          + spec_array_left[i].sum_speciation_edges_subtree
          + spec_array_right[j].sum_speciation_edges_subtree;
        int combined_spec_num = data->num_known_speciation_edges
          + i + j;
        double combined_val =
          spec_array_left[i].coalescent_sum_subtree
          + spec_array_right[j].coalescent_sum_subtree
          + compute_speciation_value(combined_spec_num, combined_spec_sum);
        if (combined_val > spec_array_act[i+j+2].score)
        {
          spec_array_act[i+j+2].score = combined_val;
          spec_array_act[i+j+2].sum_speciation_edges_subtree =
            spec_array_left[i].sum_speciation_edges_subtree
            + spec_array_right[j].sum_speciation_edges_subtree
            + tree->left->length + tree->right->length;
          spec_array_act[i+j+2].coalescent_sum_subtree =
            spec_array_left[i].coalescent_sum_subtree
            + spec_array_right[j].coalescent_sum_subtree;
        }
      }
    }
  }
}

void ptp_multi_heuristic(rtree_t * tree)
{
  init_tree_data(tree);
  init_additional_tree_data(tree);
  // TODO: Implement the heuristic

  multi_traversal(tree);
  double max = 0;
  spec_entry* spec_array = ((node_information*) tree->data)->spec_array;
  int i;
  for (i = 0; i < ((node_information*)tree->data)->num_edges_subtree; i++)
  {
    if (max < spec_array[i].score)
    {
      max = spec_array[i].score;
    }
  }
  printf("Best score found: %.6f\n", max);
}
