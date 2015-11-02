/*
    Copyright (C) 2015 Tomas Flouri

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

static int indend_space = 4;

static void print_node_info(rtree_t * tree)
{
  printf (" %s", tree->label);
  printf (" %f", tree->length);
  printf("\n");
}

static void print_tree_recurse(rtree_t * tree, 
                               int indend_level, 
                               int * active_node_order)
{
  int i,j;

  if (!tree) return;

  for (i = 0; i < indend_level; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indend_space-1; ++j)
      printf(" ");
  }
  printf("\n");

  for (i = 0; i < indend_level-1; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indend_space-1; ++j)
      printf(" ");
  }

  printf("+");
  for (j = 0; j < indend_space-1; ++j)
    printf ("-");
  if (tree->left || tree->right) printf("+");

  print_node_info(tree);

  if (active_node_order[indend_level-1] == 2) 
    active_node_order[indend_level-1] = 0;

  active_node_order[indend_level] = 1;
  print_tree_recurse(tree->left,
                     indend_level+1,
                     active_node_order);
  active_node_order[indend_level] = 2;
  print_tree_recurse(tree->right,
                     indend_level+1,
                     active_node_order);

}

static int tree_indend_level(rtree_t * tree, int indend)
{
  if (!tree) return indend;

  int a = tree_indend_level(tree->left,  indend+1);
  int b = tree_indend_level(tree->right, indend+1);

  return (a > b ? a : b);
}

void rtree_show_ascii(rtree_t * tree)
{
  
  int indend_max = tree_indend_level(tree,0);

  int * active_node_order = (int *)malloc((size_t)(indend_max+1) * sizeof(int));
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  print_node_info(tree);
  print_tree_recurse(tree->left,  1, active_node_order);
  active_node_order[0] = 2;
  print_tree_recurse(tree->right, 1, active_node_order);
  free(active_node_order);
}

static char * rtree_export_newick_recursive(rtree_t * root)
{
  char * newick;
  char * support;

  if (!root) return NULL;

  if (!(root->left) || !(root->right))
    asprintf(&newick, "%s:%f", root->label, root->length);
  else
  {
    char * subtree1 = rtree_export_newick_recursive(root->left);
    char * subtree2 = rtree_export_newick_recursive(root->right);

    if (opt_bayes_single || opt_bayes_multi)
      asprintf(&support, "%f", root->support);

    asprintf(&newick, "(%s,%s)%s:%f", subtree1,
                                      subtree2,
                                      (opt_bayes_single || opt_bayes_multi) ? support : "",
                                      //root->label ? root->label : "",
                                      root->length);
    
    if (opt_bayes_single || opt_bayes_multi)
      free(support);

    free(subtree1);
    free(subtree2);
  }

  return newick;
}

char * rtree_export_newick(rtree_t * root)
{
  char * newick;
  char * support;

  if (!root) return NULL;

  if (!(root->left) || !(root->right))
    asprintf(&newick, "%s:%f", root->label, root->length);
  else
  {
    char * subtree1 = rtree_export_newick_recursive(root->left);
    char * subtree2 = rtree_export_newick_recursive(root->right);

    if (opt_bayes_single || opt_bayes_multi)
      asprintf(&support, "%f", root->support);

    asprintf(&newick, "(%s,%s)%s:%f;", subtree1,
                                       subtree2,
                                       (opt_bayes_single || opt_bayes_multi) ? support : "",
                                       //root->label ? root->label : "",
                                       root->length);
    if (opt_bayes_single || opt_bayes_multi)
      free(support);

    free(subtree1);
    free(subtree2);
  }

  return newick;
}

static void rtree_traverse_recursive(rtree_t * node,
                                     int (*cbtrav)(rtree_t *),
                                     int * index,
                                     rtree_t ** outbuffer)
{
  if (!node->left)
  {
    if (cbtrav(node))
    {
      outbuffer[*index] = node;
      *index = *index + 1;
    }
    return;
  }
  if (!cbtrav(node))
    return;
  rtree_traverse_recursive(node->left, cbtrav, index, outbuffer);
  rtree_traverse_recursive(node->right, cbtrav, index, outbuffer);

  outbuffer[*index] = node;
  *index = *index + 1;
}

int rtree_traverse(rtree_t * root,
                   int (*cbtrav)(rtree_t *),
                   rtree_t ** outbuffer)
{
  int index = 0;

  if (!root->left) return -1;

  /* we will traverse an unrooted tree in the following way
      
           root
            /\
           /  \
        left   right

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  rtree_traverse_recursive(root, cbtrav, &index, outbuffer);
  return index;
}

static void rtree_query_tipnodes_recursive(rtree_t * node,
                                           rtree_t ** node_list,
                                           int * index)
{
  if (!node) return;

  if (!node->left)
  {
    node_list[*index] = node;
    *index = *index + 1;
    return;
  }

  rtree_query_tipnodes_recursive(node->left,  node_list, index);
  rtree_query_tipnodes_recursive(node->right, node_list, index);
}

int rtree_query_tipnodes(rtree_t * root,
                         rtree_t ** node_list)
{
  int index = 0;

  if (!root) return 0;
  if (!root->left)
  {
    node_list[index++] = root;
    return index;
  }

  rtree_query_tipnodes_recursive(root->left,  node_list, &index);
  rtree_query_tipnodes_recursive(root->right, node_list, &index);

  return index;
}

static void rtree_query_innernodes_recursive(rtree_t * root,
                                             rtree_t ** node_list,
                                             int * index)
{
  if (!root) return;
  if (!root->left) return;

  /* postorder traversal */

  rtree_query_innernodes_recursive(root->left,  node_list, index);
  rtree_query_innernodes_recursive(root->right, node_list, index);

  node_list[*index] = root;
  *index = *index + 1;
  return;
}

int rtree_query_innernodes(rtree_t * root,
                           rtree_t ** node_list)
{
  int index = 0;

  if (!root) return 0;
  if (!root->left) return 0;

  rtree_query_innernodes_recursive(root->left,  node_list, &index);
  rtree_query_innernodes_recursive(root->right, node_list, &index);

  node_list[index++] = root;

  return index;
}

void rtree_reset_info(rtree_t * root)
{
  if (!root->left)
  {
    root->leaves = 1;
    root->edge_count = 0;
    root->edgelen_sum = 0;
    return;
  }
  
  rtree_reset_info(root->left);
  rtree_reset_info(root->right);

  root->leaves = root->left->leaves + root->right->leaves;
  root->edge_count = root->left->edge_count + 
                     root->right->edge_count;
  root->edgelen_sum = root->left->edgelen_sum +
                      root->right->edgelen_sum;

  if (root->left->length > opt_minbr)
  {
    root->edge_count++;
    root->edgelen_sum += root->left->length;
  }
  if (root->right->length > opt_minbr)
  {
    root->edge_count++;
    root->edgelen_sum += root->right->length;
  }
}

void rtree_print_tips(rtree_t * node, FILE * out)
{
  if (node->left)  rtree_print_tips(node->left,out);
  if (node->right) rtree_print_tips(node->right,out);

  if (!node->left && !node->right)
    fprintf(out, "%s\n", node->label);
}
