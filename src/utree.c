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

static void print_node_info(utree_t * tree)
{
  printf (" %s", tree->label);
  printf (" %f", tree->length);
  printf("\n");
}

static void print_tree_recurse(utree_t * tree, 
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
  if (tree->next) printf("+");

  print_node_info(tree);

  if (active_node_order[indend_level-1] == 2) 
    active_node_order[indend_level-1] = 0;

  if (tree->next)
  {
    active_node_order[indend_level] = 1;
    print_tree_recurse(tree->next->back,
                       indend_level+1,
                       active_node_order);
    active_node_order[indend_level] = 2;
    print_tree_recurse(tree->next->next->back, 
                       indend_level+1,
                       active_node_order);
  }

}

static int tree_indend_level(utree_t * tree, int indend)
{
  if (!tree->next) return indend+1;

  int a = tree_indend_level(tree->next->back,       indend+1);
  int b = tree_indend_level(tree->next->next->back, indend+1);

  return (a > b ? a : b);
}

void utree_show_ascii(utree_t * tree)
{
  int a, b;
  
  a = tree_indend_level(tree->back,1);
  b = tree_indend_level(tree,0);
  int max_indend_level = (a > b ? a : b);


  int * active_node_order = (int *)malloc((max_indend_level+1) * sizeof(int));
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  print_tree_recurse(tree->back,             1, active_node_order);
  print_tree_recurse(tree->next->back,       1, active_node_order);
  active_node_order[0] = 2;
  print_tree_recurse(tree->next->next->back, 1, active_node_order);
  free(active_node_order);
}

static char * newick_utree_recurse(utree_t * root)
{
  char * newick;

  if (!root->next)
    asprintf(&newick, "%s:%f", root->label, root->length);
  else
  {
    char * subtree1 = newick_utree_recurse(root->next->back);
    char * subtree2 = newick_utree_recurse(root->next->next->back);

    asprintf(&newick, "(%s,%s)%s:%f", subtree1,
                                      subtree2,
                                      root->label ? root->label : "",
                                      root->length);
    free(subtree1);
    free(subtree2);
  }

  return newick;
}

char * utree_export_newick(utree_t * root)
{
  char * newick;

  if (!root) return NULL;

  char * subtree1 = newick_utree_recurse(root->back);
  char * subtree2 = newick_utree_recurse(root->next->back);
  char * subtree3 = newick_utree_recurse(root->next->next->back);

  asprintf(&newick, "(%s,%s,%s)%s:%f;", subtree1,
                                        subtree2,
                                        subtree3,
                                        root->label ? root->label : "",
                                        root->length);
  free(subtree1);
  free(subtree2);
  free(subtree3);

  return (newick);

}

static void utree_traverse_recursive(utree_t * node,
                                     int (*cbtrav)(utree_t *),
                                     int * index,
                                     utree_t ** outbuffer)
{
  if (!node->next)
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

  utree_traverse_recursive(node->next->back, cbtrav, index, outbuffer);
  utree_traverse_recursive(node->next->next->back, cbtrav, index, outbuffer);

  outbuffer[*index] = node;
  *index = *index + 1;
}

int utree_traverse(utree_t * root,
                   int (*cbtrav)(utree_t *),
                   utree_t ** outbuffer)
{
  int index = 0;

  if (!root->next) return -1;

  /* we will traverse an unrooted tree in the following way
      
              2
            /
      1  --*
            \
              3

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  utree_traverse_recursive(root->back, cbtrav, &index, outbuffer);
  utree_traverse_recursive(root, cbtrav, &index, outbuffer);

  return index;
}


static void utree_query_tipnodes_recursive(utree_t * node,
                                           utree_t ** node_list,
                                           int * index)
{
  if (!node->next)
  {
    node_list[*index] = node;
    *index = *index + 1;
    return;
  }

  utree_query_tipnodes_recursive(node->next->back, node_list, index);
  utree_query_tipnodes_recursive(node->next->next->back, node_list, index);
}

int utree_query_tipnodes(utree_t * root,
                         utree_t ** node_list)
{
  int index = 0;

  if (!root) return 0;

  if (!root->next) root = root->back;

  utree_query_tipnodes_recursive(root->back, node_list, &index);

  utree_query_tipnodes_recursive(root->next->back, node_list, &index);
  utree_query_tipnodes_recursive(root->next->next->back, node_list, &index);

  return index;
}

static void utree_query_innernodes_recursive(utree_t * node,
                                             utree_t ** node_list,
                                             int * index)
{
  if (!node->next) return;

  /* postorder traversal */

  utree_query_innernodes_recursive(node->next->back, node_list, index);
  utree_query_innernodes_recursive(node->next->next->back, node_list, index);

  node_list[*index] = node;
  *index = *index + 1;
  return;
}

int utree_query_innernodes(utree_t * root,
                           utree_t ** node_list)
{
  int index = 0;

  if (!root) return 0;
  if (!root->next) root = root->back;

  utree_query_innernodes_recursive(root->back, node_list, &index);

  utree_query_innernodes_recursive(root->next->back, node_list, &index);
  utree_query_innernodes_recursive(root->next->next->back, node_list, &index);

  node_list[index++] = root;

  return index;
}

static rtree_t * utree_rtree(utree_t * unode)
{
  rtree_t * rnode = (rtree_t *)xmalloc(sizeof(rtree_t)); 

  rnode->event = EVENT_COALESCENT;

  if (unode->label)
    rnode->label = strdup(unode->label);
  else
    rnode->label = NULL;
  rnode->length = unode->length;

  if (!unode->next) 
  {
    rnode->left = NULL;
    rnode->right = NULL;
    return rnode;
  }

  rnode->left = utree_rtree(unode->next->back);
  rnode->right = utree_rtree(unode->next->next->back);

  rnode->left->parent = rnode;
  rnode->right->parent = rnode;

  return rnode;
}

static utree_t * find_longest_branchtip(utree_t ** node_list, int tip_count)
{
  int i;
  
  double branch_length = 0;
  int index = 0;

  /* check whether there exists a tip with the outgroup label */
  for (i = 0; i < tip_count; ++i)
    if (node_list[i]->length > branch_length)
    {
      index = i;
      branch_length = node_list[i]->length;
    }

  fprintf(stdout, 
          "Selected %s as outgroup based on longest tip-branch criterion\n",
          node_list[index]->label);

  return node_list[index];
}

static utree_t * find_outgroup_node(utree_t ** node_list, int tip_count)
{
  int i;

  /* check whether there exists a tip with the outgroup label */
  for (i = 0; i < tip_count; ++i)
    if (!strcmp(node_list[i]->label, opt_outgroup)) break;

  if (i == tip_count)
    fatal("Outgroup not among tips");

  fprintf(stdout, "Rooting with outgroup %s\n", node_list[i]->label);

  return node_list[i];
}

static utree_t * find_outgroup_mrca(utree_t ** node_list,
                                    utree_t * root,
                                    int tip_count)
{
  int i,j,k,count;
  utree_t * outgroup;

  size_t taxon1_len = strcspn(opt_outgroup, ",");

  char * taxon1 = strndup(opt_outgroup, taxon1_len);
  char * taxon2 = strdup(opt_outgroup+taxon1_len+1);

  fprintf(stdout, "Rooting with outgroups %s and %s\n", taxon1, taxon2);

  j = k = tip_count;
  for (count = 0, i = 0; i < tip_count; ++i)
  {
    if (!strcmp(node_list[i]->label, taxon1))
    {
      count++;
      j = i;
    }
    if (!strcmp(node_list[i]->label, taxon2))
    {
      k = i;
      count++;
    }
    if (count == 2) break;
  }

  if (j == tip_count)
    fatal("Outgroup %s not among tips", taxon1);

  if (k == tip_count)
    fatal("Outgroup %s not among tips", taxon2);

  free(taxon1);
  free(taxon2);

  /* find LCA */
  lca_init(root);
  outgroup = lca_compute(node_list[j], node_list[k]);
  lca_destroy();

  return outgroup;
}

rtree_t * utree_convert_rtree(utree_t * root, int tip_count)
{
  utree_t * outgroup;

  /* query tip nodes */
  utree_t ** node_list = (utree_t **)xmalloc(tip_count * sizeof(utree_t *));
  utree_query_tipnodes(root, node_list);

  /* find outgroup */
  if (!opt_outgroup)
    outgroup = find_longest_branchtip(node_list, tip_count);
  else if (!strchr(opt_outgroup, ','))
    outgroup = find_outgroup_node(node_list, tip_count);
  else
    outgroup = find_outgroup_mrca(node_list, root, tip_count);

  rtree_t * rnode = (rtree_t *)xmalloc(sizeof(rtree_t));
  rnode->left = utree_rtree(outgroup);
  rnode->right = utree_rtree(outgroup->back);
  rnode->parent = NULL;

  rnode->left->parent = rnode;
  rnode->right->parent = rnode;
  rnode->left->length /= 2;
  rnode->right->length /= 2;
  rnode->label = NULL;
  rnode->length = 0;

  /* reset per-node leaves and valid edges */
  rtree_reset_info(rnode);

  free(node_list);

  return rnode;
}
