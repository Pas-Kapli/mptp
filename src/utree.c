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


  int * active_node_order = (int *)malloc((size_t)(max_indend_level+1) * 
                                          sizeof(int));
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

static void utree_traverse_postorder_recursive(utree_t * node,
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

  utree_traverse_postorder_recursive(node->next->back, cbtrav, index, outbuffer);
  utree_traverse_postorder_recursive(node->next->next->back, cbtrav, index, outbuffer);

  if (cbtrav(node))
  {
    outbuffer[*index] = node;
    *index = *index + 1;
  }
}

static int cb_outgroup(utree_t * node)
{
  /* if it's a tip */
  if (!node->next)
    return 0;

  /* if inner node */
  if (node->next->back->mark == 1 || node->next->next->back->mark == 1)
    node->mark = 1;
  else
    node->mark = 0;

  node->next->mark = node->next->back->mark;
  node->next->next->mark = node->next->next->back->mark;

  return node->mark;
}

int utree_traverse_postorder(utree_t * root,
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

  utree_traverse_postorder_recursive(root->back, cbtrav, &index, outbuffer);
  utree_traverse_postorder_recursive(root, cbtrav, &index, outbuffer);

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
    rnode->label = xstrdup(unode->label);
  else
    rnode->label = NULL;
  rnode->length = unode->length;
  rnode->data = NULL;
  rnode->mark = 0;

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

utree_t * utree_longest_branchtip(utree_t * node, int tip_count)
{
  int i;
  int index = 0;
  double branch_length = 0;
  utree_t * outgroup = NULL;

  /* query tip nodes */
  utree_t ** tip_nodes_list = (utree_t **)xmalloc((size_t)tip_count * sizeof(utree_t *));
  utree_query_tipnodes(node, tip_nodes_list);

  for (i = 0; i < tip_count; ++i)
    if (tip_nodes_list[i]->length > branch_length)
    {
      index = i;
      branch_length = tip_nodes_list[i]->length;
    }
  
  outgroup = tip_nodes_list[index];

  free(tip_nodes_list);

  return outgroup;
}

rtree_t * utree_crop(utree_t * lca)
{
  /* is the back of the lca a tip? */
  if (!lca->back->next)
    return NULL; 

  rtree_t * root = (rtree_t *)xmalloc(sizeof(rtree_t));
  
  /* clone the two subtrees */
  root->left  = utree_rtree(lca->back->next->back);
  root->right = utree_rtree(lca->back->next->next->back);

  root->parent = NULL;
  root->length = 0;
  root->label  = NULL;
  root->data   = NULL;
  root->mark   = 0;

  root->left->parent = root;
  root->right->parent = root;

  rtree_reset_info(root);

  return root;
}

rtree_t * utree_convert_rtree(utree_t * outgroup)
{
  rtree_t * root = (rtree_t *)xmalloc(sizeof(rtree_t));
  root->left   = utree_rtree(outgroup);
  root->right  = utree_rtree(outgroup->back);

  root->left->parent  = root;
  root->right->parent = root;

  root->left->length /= 2;
  root->right->length /= 2;

  root->label  = NULL;
  root->length = 0;
  root->parent = NULL;
  root->event  = EVENT_COALESCENT;
  root->data   = NULL;
  root->mark   = 0;

  /* reset per-node leaves and valid edges */
  rtree_reset_info(root);

  return root;
  
}

static utree_t ** utree_tipstring_nodes(utree_t * root,
                                        char * tipstring,
                                        int utree_tip_count,
                                        unsigned int * tiplist_count)
{
  unsigned int i;
  unsigned int k;
  unsigned int commas_count = 0;

  char * taxon;
  unsigned int taxon_len;

  ENTRY * found = NULL;

  for (i = 0; i < strlen(tipstring); ++i)
    if (tipstring[i] == ',')
      commas_count++;
  
  utree_t ** node_list = (utree_t **)xmalloc(utree_tip_count*sizeof(utree_t *));
  utree_query_tipnodes(root, node_list);

  utree_t ** out_node_list = (utree_t **)xmalloc((commas_count+1) *
                                                   sizeof(utree_t *));

  /* create a hashtable of tip labels */
  hcreate(2 * utree_tip_count);

  for (i = 0; i < (unsigned int)utree_tip_count; ++i)
  {
    ENTRY entry;
    entry.key  = node_list[i]->label;
    entry.data = node_list[i];
    hsearch(entry,ENTER);
  }

  char * s = tipstring;
  
  k = 0;
  while (*s)
  {
    /* get next tip */
    taxon_len = strcspn(s, ",");
    if (!taxon_len)
      fatal("Erroneous prune list format (double comma)/taxon missing");

    taxon = strndup(s, taxon_len);

    /* search tip in hash table */
    ENTRY query;
    query.key = taxon;
    found = NULL;
    found = hsearch(query,FIND);
    
    if (!found)
      fatal("Taxon %s does not appear in the tree", taxon);

    /* store pointer in output list */
    out_node_list[k++] = (utree_t *)(found->data);

    /* free tip label, and move to the beginning of next tip if available */
    free(taxon);
    s += taxon_len;
    if (*s == ',') 
      s += 1;
  }

  /* kill the hash table */
  hdestroy();

  free(node_list);

  /* return number of tips in the list */
  *tiplist_count = commas_count + 1;

  /* return tip node list */
  return out_node_list;
}

utree_t * utree_lca(utree_t * root,
                    utree_t ** tip_nodes,
                    unsigned int count,
                    unsigned int utree_tip_count)
{
  long i;
  utree_t * lca;
  utree_t ** path;

  /* allocate a path */
  path = (utree_t **)xmalloc((size_t)utree_tip_count *
                                  sizeof(utree_t **));

  /* mark all tip nodes */
  for (i = 0; i < count; ++i)
    tip_nodes[i]->mark = 1;

  /* traverse the tree with the cb_outgroup callback to get the inner nodes
     of the subtree formed by the outgroup */
  int path_len = utree_traverse_postorder(tip_nodes[0]->back,
                                          cb_outgroup,
                                          path);


  /* there must be exactly one inner node that does not have all three 
     directions mark. That one will be the root of the outgroup subtree */
  int root_count = 0;
  for (i = 0; i < path_len; ++i)
    if (!(path[i]->mark && path[i]->next->mark && path[i]->next->next->mark))
    {
      root_count++;
      lca = path[i];
    }

  /* deallocate path */
  free(path);

  /* if we had more than one inner nodes with less than three directions marked
     then not all tips of a subtree were specified (invalid outgroup) */
  if (root_count != 1) return NULL;
  while (lca->mark == 1) lca = lca->next;

  /* return the LCA */
  return lca;
}

utree_t * utree_outgroup_lca(utree_t * root, unsigned int tip_count)
{
  unsigned int og_tips_count;
  utree_t * og_root;
  utree_t ** og_tips;

  /* get all nodes that have labels equal to the comma separated string in
     opt_outgroup */
  og_tips = utree_tipstring_nodes(root,
                                  opt_outgroup,
                                  tip_count,
                                  &og_tips_count);

  if (og_tips_count == 1)
  {
    og_root = og_tips[0];
  }
  else
  {
    /* find the LCA of the tips in og_tips. Note that, *all* tips of the desired
       subtree *must* be specified */
    og_root = utree_lca(root, og_tips, og_tips_count, tip_count);
  }

  free(og_tips);

  /* return the LCA (root of the outgroup subtree */
  return og_root;
}
