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
  info->is_real_mrca = false;
  info->is_input_mrca = false;
  info->current_species_real = -1;
  info->current_species_input = -1;
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
  tree->data = NULL;
}

void identify_alternative_taxa(char * scorefile, int num_leaves,
    rtree_t ** leaves)
{
  FILE * scorefile_in = fopen(scorefile, "r");
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  if (!scorefile_in)
  {
    snprintf(errmsg, 200, "Unable to open file (%s)", scorefile);
  }
  else
  {
    bool collecting = false;
    int current_species_idx = 0;
    while ((read = getline(&line, &len, scorefile_in)) != -1)
    {
      //printf("line: %s", line);
      if (strstr(line, "Species")!= NULL)
      {
        //printf("A new species starts.\n");
        collecting = true;
        current_species_idx++;
      }
      else if (collecting &&
        ((strcmp(line,"\n")==0
      || strstr(line, "Writing tree file") != NULL)
        || strstr(line, "Number of") != NULL))
      {
        //printf("A species has ended.\n");
        collecting = false;
      }
      else if (collecting)
      {
        // find leaf node with the same label as line,
        // give it current_species_input = current_species_idx

        bool found = false;
        int i;
        for (i = 0; i < num_leaves; i++)
        {
          strtok(line,"\n");
          if (strcmp(leaves[i]->label, line) == 0)
          {
            ((score_information*) (leaves[i]->data))->current_species_input
                = current_species_idx;
            found = true;
            break;
          }
        }
        if (!found)
        {
          printf("Not found: %s", line);
        }
      }
    }
    fclose(scorefile_in);
    if (line)
    {
      free(line);
    }
  }
}

void retrieve_mrca_nodes(rtree_t * tree, int * num_species_real,
  int * num_species_input)
{
  // Do a post-order traversal for finding common ancestor nodes.
  // As long as the current letters are the same, it is the same species.
  if (tree->left)
  {
    retrieve_mrca_nodes(tree->left, num_species_real, num_species_input);
  }
  if (tree->right)
  {
    retrieve_mrca_nodes(tree->right, num_species_real, num_species_input);
  }
  if (tree->left && tree->right) // inner node
  {
    score_information * data_current = (score_information*) tree->data;
    score_information * data_left = (score_information*) tree->left->data;
    score_information * data_right = (score_information*) tree->right->data;

    if (data_left->current_species_real == data_right->current_species_real)
    {
      data_current->current_species_real = data_left->current_species_real;
    }
    else
    {
      if (data_left->current_species_real != -1)
      {
        data_left->marked = true;
        data_left->is_real_mrca = true;
        (*num_species_real)++;
      }
      if (data_right->current_species_real != -1)
      {
        data_right->marked = true;
        data_right->is_real_mrca = true;
        (*num_species_real)++;
      }
    }

    if (data_left->current_species_input == data_right->current_species_input)
    {
      data_current->current_species_input = data_left->current_species_input;
    }
    else
    {
      if (data_left->current_species_input != -1)
      {
        data_left->marked = true;
        data_left->is_input_mrca = true;
        (*num_species_input)++;
      }
      if (data_right->current_species_input != -1)
      {
        data_right->marked = true;
        data_right->is_input_mrca = true;
        (*num_species_input)++;
      }
    }
  }
  else if (!tree->left && !tree->right) // leaf node
  {
    char * label = tree->label;
    char * species_idx_real = strtok(label,".");
    ((score_information*) tree->data)->current_species_real =
        atoi(species_idx_real);
  }
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

void compute_tree_penalty_score(rtree_t * current_node, rtree_t * root,
  int * score_ptr)
{
  // Go from each marked node up towards the root, if another marked node is
  // encountered, add up the number of steps taken to get there to the penalty.
  // If no other marked node was in the way, add nothing to the penalty.

  if (current_node->left)
  {
    compute_tree_penalty_score(current_node->left, root, score_ptr);
  }
  if (current_node->right)
  {
    compute_tree_penalty_score(current_node->right, root, score_ptr);
  }

  if (((score_information*)current_node->data)->marked)
  { // current node is a mrca node
    *(score_ptr) += walk_to_root(current_node, root);
  }
}

void collect_mrca_nodes_recursive(rtree_t * tree, rtree_t ** mrca_real_list,
  rtree_t ** mrca_input_list, int * index_real, int * index_input)
{
  if (tree->left)
  {
    collect_mrca_nodes_recursive(tree->left, mrca_real_list, mrca_input_list,
      index_real, index_input);
  }
  if (tree->right)
  {
    collect_mrca_nodes_recursive(tree->right, mrca_real_list, mrca_input_list,
      index_real, index_input);
  }

  score_information * data_current = (score_information*) tree->data;
  if (data_current->is_real_mrca)
  {
    mrca_real_list[(*index_real)] = tree;
    (*index_real)++;
  }
  if (data_current->is_input_mrca)
  {
    mrca_input_list[(*index_input)] = tree;
    (*index_input)++;
  }
}

void collect_mrca_nodes(rtree_t * tree, rtree_t ** mrca_real_list,
  rtree_t ** mrca_input_list)
{
  int index_real = 0;
  int index_input = 0;
  collect_mrca_nodes_recursive(tree, mrca_real_list, mrca_input_list,
    &index_real, &index_input);
}

double compute_entropy(rtree_t ** mrca_list, int num_species, int num_taxa)
{
  double entropy = 0;
  int i;
  for (i = 0; i < num_species; i++)
  {
    double percentage = (double) mrca_list[i]->leaves / (double) num_taxa;
    entropy += percentage * log(percentage);
  }
  entropy *= -1;
  return entropy;
}

int count_common_taxa(rtree_t * mrca_one, rtree_t * mrca_two)
{
  // get leaves
  rtree_t ** leaves_one_list = calloc(mrca_one->leaves, sizeof(rtree_t));
  rtree_t ** leaves_two_list = calloc(mrca_two->leaves, sizeof(rtree_t));
  rtree_query_tipnodes(mrca_one, leaves_one_list);
  rtree_query_tipnodes(mrca_two, leaves_two_list);

  int num_common_taxa = 0;
  // TODO: Make this search more efficient.
  int i;
  int j;
  for (i = 0; i < mrca_one->leaves; i++)
  {
    for (j = 0; j < mrca_two->leaves; j++)
    {
      if (leaves_one_list[i] == leaves_two_list[j])
      {
        num_common_taxa++;
      }
    }
  }
  num_common_taxa /= 2;

  free(leaves_one_list);
  free(leaves_two_list);
  return num_common_taxa;
}

double compute_nmi_score(rtree_t ** mrca_real_list, int num_species_real,
  rtree_t ** mrca_input_list, int num_species_input, int num_taxa)
{
  double entropy_real = compute_entropy(mrca_real_list, num_species_real,
    num_taxa);
  double entropy_input = compute_entropy(mrca_input_list, num_species_input,
    num_taxa);
  double max_entropy = entropy_real;
  if (entropy_input > entropy_real)
  {
    max_entropy = entropy_input;
  }

  // compute mutual information
  double mutual_information = 0;
  int i;
  int j;
  for (i = 0; i < num_species_real; i++)
  {
    for (j = 0; j < num_species_input; j++)
    {
      rtree_t * mrca_real = mrca_real_list[i];
      rtree_t * mrca_input = mrca_input_list[j];
      double num_taxa_real = (double) mrca_real->leaves;
      double num_taxa_input = (double) mrca_input->leaves;
      double num_common_taxa = (double) count_common_taxa(mrca_real, mrca_input);

      /*printf("num_taxa: %d\n", num_taxa);
      printf("num_taxa_real: %d\n", num_taxa_real);
      printf("num_taxa_input: %d\n", num_taxa_input);
      printf("num_common_taxa: %d\n", num_common_taxa);*/

      if (num_common_taxa != 0)
      {
        mutual_information += (num_common_taxa / num_taxa)
          * log((num_common_taxa / num_taxa)
                 / (num_taxa_real * num_taxa_input / (num_taxa * num_taxa)));
      }
    }
  }

  printf("Mutual information: %.6f\n", mutual_information);
  printf("Entropy real: %.6f\n", entropy_real);
  printf("Entropy input: %.6f\n", entropy_input);
  printf("Maximum Entropy: %.6f\n", max_entropy);
  return 1 - (mutual_information / max_entropy);
}

void compute_score(rtree_t * tree, rtree_t ** mrca_list, int num_species,
  double * score_single, double * score_multi)
{
  node_information * data_root = (node_information*) (tree->data);
  double speciation = 0;
  double coalescent_single = 0;
  double coalescent_multi = 0;
  int num_edges_coalescent = 0;
  double sum_edges_coalescent = 0;
  int i;
  for (i = 0; i < num_species; i++)
  {
    node_information * data_current = (node_information*) (mrca_list[i]->data);
    coalescent_multi += data_current->coalescent;
    //printf("  Coalescent multi is now adding: %.6f\n", data_current->coalescent);
    num_edges_coalescent += data_current->num_edges_subtree;
    sum_edges_coalescent += data_current->sum_edges_subtree;
  }
  speciation = compute_loglikelihood(
    data_root->num_edges_subtree - num_edges_coalescent,
    data_root->sum_edges_subtree - sum_edges_coalescent);
  coalescent_single = compute_loglikelihood(num_edges_coalescent,
    sum_edges_coalescent);

  (*score_single) = coalescent_single + speciation;
  (*score_multi) = coalescent_multi + speciation;
}

void score_delimitation_tree(char * scorefile, rtree_t * tree, double min_br)
{
  rtree_t ** leaves_list = calloc(tree->leaves, sizeof(rtree_t));
  rtree_query_tipnodes(tree, leaves_list);

  init_tree_data_score(tree);
  // Identify taxa of the scorefile solution
  identify_alternative_taxa(scorefile, tree->leaves, leaves_list);

  // Mark mrca nodes
  int num_species_real = 0;
  int num_species_input = 0;
  retrieve_mrca_nodes(tree, &num_species_real, &num_species_input);
  printf("Number of real species: %d\n", num_species_real);
  printf("Number of species in input file: %d\n", num_species_input);

  rtree_t ** mrca_real_list = calloc(num_species_real, sizeof(rtree_t));
  rtree_t ** mrca_input_list = calloc(num_species_input, sizeof(rtree_t));
  collect_mrca_nodes(tree, mrca_real_list, mrca_input_list);

  int tree_penalty_score = 0;
  compute_tree_penalty_score(tree, tree, &tree_penalty_score);
  printf("Tree penalty score: %d\n", tree_penalty_score);
  double nmi_score = compute_nmi_score(mrca_real_list, num_species_real,
    mrca_input_list, num_species_input, tree->leaves);
  printf("NMI score: %.6f\n", nmi_score);

  free_tree_data_score(tree);

  init_tree_data(tree, min_br);
  double score_single_input = 0;
  double score_multi_input = 0;
  double score_single_real = 0;
  double score_multi_real = 0;
  compute_score(tree, mrca_input_list, num_species_input,
    &score_single_input, &score_multi_input);
  compute_score(tree, mrca_real_list, num_species_real,
    &score_single_real, &score_multi_real);
  printf("Score input single: %.6f\n", score_single_input);
  printf("Score input multi: %.6f\n", score_multi_input);
  printf("Score real single: %.6f\n", score_single_real);
  printf("Score real multi: %.6f\n", score_multi_real);
  free_tree_data(tree);

  free(leaves_list);
  free(mrca_real_list);
  free(mrca_input_list);
}
