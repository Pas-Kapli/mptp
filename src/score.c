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
#include <string.h>

void init_tree_data_score(rtree_t * tree)
{
  if (tree->left)
  {
    init_tree_data_score(tree->left);
  }
  if (tree->right)
  {
    init_tree_data_score(tree->right);
  }
  score_information * info = malloc(sizeof(score_information));
  info->marked = false;
  info->current_species = -1;
  tree->data = info;
}

void free_tree_data_score(rtree_t * tree)
{
  if (tree->left)
  {
    free_tree_data_score(tree->left);
  }
  if (tree->right)
  {
    free_tree_data_score(tree->right);
  }
  free((score_information*) (tree->data));
}

void retrieve_optimal_mrca_nodes(rtree_t * tree)
{
  // Do a post-order traversal for finding common ancestor nodes.
  // As long as the current letters are the same, it is the same species.
  if (tree->left)
  {
    retrieve_optimal_mrca_nodes(tree->left);
  }
  if (tree->right)
  {
    retrieve_optimal_mrca_nodes(tree->right);
  }
  if (tree->left && tree->right) // inner node
  {
    score_information * data_current = (score_information*) tree->data;
    score_information * data_left = (score_information*) tree->left->data;
    score_information * data_right = (score_information*) tree->right->data;

    if (data_left->current_species == data_right->current_species)
    {
      data_current->current_species = data_left->current_species;
    }
    else
    {
      if (data_left->current_species != -1)
      {
        data_left->marked = true;
        printf("Found mrca for species %d\n", data_left->current_species);
      }
      if (data_right->current_species != -1)
      {
        data_right->marked = true;
        printf("Found mrca for species %d\n", data_right->current_species);
      }
    }
  }
  else if (!tree->left && !tree->right) // leaf node
  {
    char * label = tree->label;
    char * species_idx = strtok(label,".");
    ((score_information*) tree->data)->current_species = atoi(species_idx);
  }
}

void retrieve_alternative_mrca_nodes(char * scorefile, rtree_t * tree)
{
  // TODO: implement this function
}

int walk_to_root(rtree_t * starting_node, rtree_t * root)
{
  int steps = 0;
  int penalty = 0;
  rtree_t * current_node = starting_node;
  while(current_node != root)
  {
    current_node = current_node->parent;
    steps++;
    if (((score_information*)current_node->data)->marked)
    {
      penalty = steps;
      break;
    }
  }
  return penalty;
}

void compute_score(rtree_t * current_node, rtree_t * root, int * score_ptr)
{
  // Go from each marked node up towards the root, if another marked node is
  // encountered, add up the number of steps taken to get there to the penalty.
  // If no other marked node was in the way, add nothing to the penalty.

  if (current_node->left)
  {
    compute_score(current_node->left, root, score_ptr);
  }
  if (current_node->right)
  {
    compute_score(current_node->right, root, score_ptr);
  }

  if (((score_information*)current_node->data)->marked) // current node is a mrca node
  {
    *(score_ptr) += walk_to_root(current_node, root);
  }
}

void score_delimitation_tree(char * scorefile, rtree_t * tree)
{
  init_tree_data_score(tree);
  // Mark mrca nodes of the optimal solution
  retrieve_optimal_mrca_nodes(tree);
  // Mark mrca nodes of the scorefile solution
  retrieve_alternative_mrca_nodes(scorefile, tree);

  int score = 0;
  compute_score(tree, tree, &score);
  printf("Tree penalty score: %d\n", score);
  free_tree_data_score(tree);
}

