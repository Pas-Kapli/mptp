/*
    Copyright (C) 2015-2017 Tomas Flouri

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

void rnode_destroy(rnode_t * root)
{
  if (!root) return;

  rnode_destroy(root->left);
  rnode_destroy(root->right);
  if (root->data)
    free(root->data);

  free(root->label);
  free(root);
}

void rtree_destroy(rtree_t * tree)
{
  if (!tree) return;
  rnode_destroy(tree->root);
  free(tree->nodes);
  free(tree);
}

rtree_t * rtree_clone(rtree_t * tree)
{
  rtree_t * clone = (rtree_t *)xcalloc(1, sizeof(rtree_t));
  clone->root = rnode_clone(tree->root, NULL);
  clone->tip_count = tree->tip_count;
  clone->inner_count = tree->inner_count;
  clone->edge_count = tree->edge_count;
  clone->nodes = (rnode_t **)xmalloc(
      (size_t)(clone->tip_count + clone->inner_count) * sizeof(rnode_t *));
  rnode_query_tipnodes(clone->root, clone->nodes);
  rnode_query_innernodes(clone->root, clone->nodes + clone->tip_count);
  return clone;
}

static int indend_space = 4;

static void print_node_info(rnode_t * tree)
{
  printf (" %s", tree->label);
  printf (" %f", tree->length);
  printf("\n");
}

static void print_tree_recurse(rnode_t * tree,
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

static int tree_indend_level(rnode_t * tree, int indend)
{
  if (!tree) return indend;

  int a = tree_indend_level(tree->left,  indend+1);
  int b = tree_indend_level(tree->right, indend+1);

  return (a > b ? a : b);
}

void rtree_show_ascii(rtree_t * tree)
{
  rnode_t * root = tree->root;

  int indend_max = tree_indend_level(root,0);

  int * active_node_order = (int *)malloc((size_t)(indend_max+1) * sizeof(int));
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  print_node_info(root);
  print_tree_recurse(root->left,  1, active_node_order);
  active_node_order[0] = 2;
  print_tree_recurse(root->right, 1, active_node_order);
  free(active_node_order);
}

static char * rnode_export_newick_recursive(rnode_t * root)
{
  char * newick;
  char * support = NULL;

  if (!root) return NULL;

  if (!(root->left) || !(root->right))
  {
    if (asprintf(&newick, "%s:%f", root->label, root->length) == -1)
      fatal("Unable to allocate enough memory.");
  }
  else
  {
    char * subtree1 = rnode_export_newick_recursive(root->left);
    char * subtree2 = rnode_export_newick_recursive(root->right);

    if (opt_mcmc)
      if (asprintf(&support, "%f", root->support) == -1)
        fatal("Unable to allocate enough memory.");

    if (asprintf(&newick, "(%s,%s)%s:%f", subtree1,
                                      subtree2,
                                      (opt_mcmc) ? support : "",
                                      root->length) == -1)
      fatal("Unable to allocate enough memory.");

    if (opt_mcmc)
      free(support);

    free(subtree1);
    free(subtree2);
  }

  return newick;
}

char * rtree_export_newick(rtree_t * tree)
{
  rnode_t * root = tree->root;
  char * newick;
  char * support = NULL;

  if (!root) return NULL;

  if (!(root->left) || !(root->right))
  {
    if (asprintf(&newick, "%s:%f", root->label, root->length) == -1)
      fatal("Unable to allocate enough memory.");
  }
  else
  {
    char * subtree1 = rnode_export_newick_recursive(root->left);
    char * subtree2 = rnode_export_newick_recursive(root->right);

    if (opt_mcmc)
      if (asprintf(&support, "%f", root->support) == -1)
        fatal("Unable to allocate enough memory.");

    if (asprintf(&newick, "(%s,%s)%s:%f;", subtree1,
                                       subtree2,
                                       (opt_mcmc) ? support : "",
                                       root->length) == -1)
      fatal("Unable to allocate enough memory.");
    if (opt_mcmc)
      free(support);

    free(subtree1);
    free(subtree2);
  }

  return newick;
}

static void rnode_traverse_recursive(rnode_t * node,
                                     int (*cbtrav)(rnode_t *),
                                     int * index,
                                     unsigned short * rstate,
                                     rnode_t ** outbuffer)
{
  double rand_double = 0;

  if (!node->left)
  {
    if (!cbtrav(node))
    {
      outbuffer[*index] = node;
      *index = *index + 1;
    }
    return;
  }
  if (!cbtrav(node))
  {
    outbuffer[*index] = node;
    *index = *index + 1;
    return;
  }

  rand_double = mptp_erand48(rstate);
  if (rand_double >= 0.5)
  {
    rnode_traverse_recursive(node->left, cbtrav, index, rstate, outbuffer);
    rnode_traverse_recursive(node->right, cbtrav, index, rstate, outbuffer);
  }
  else
  {
    rnode_traverse_recursive(node->right, cbtrav, index, rstate, outbuffer);
    rnode_traverse_recursive(node->left, cbtrav, index, rstate, outbuffer);
  }

}

int rnode_traverse(rnode_t * root,
                   int (*cbtrav)(rnode_t *),
                   unsigned short * rstate,
                   rnode_t ** outbuffer)
{
  int index = 0;

  if (!root->left) return -1;

  /* we will traverse an rooted tree in the following way

           root
            /\
           /  \
        left   right

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  rnode_traverse_recursive(root, cbtrav, &index, rstate, outbuffer);
  return index;
}

static void rnode_traverse_postorder_recursive(rnode_t * node,
                                               int (*cbtrav)(rnode_t *),
                                               int * index,
                                               rnode_t ** outbuffer)
{
  if (!node) return;

  rnode_traverse_postorder_recursive(node->left,  cbtrav, index, outbuffer);
  rnode_traverse_postorder_recursive(node->right, cbtrav, index, outbuffer);

  if (cbtrav(node))
  {
    outbuffer[*index] = node;
    *index = *index + 1;
  }
}


int rnode_traverse_postorder(rnode_t * root,
                             int (*cbtrav)(rnode_t *),
                             rnode_t ** outbuffer)
{
  int index = 0;

  if (!root->left) return -1;

  /* we will traverse an unrooted tree in the following way

           root
            /\
           /  \
        left   right

     at each node the callback function is called to decide whether to
     place the node in the list */

  rnode_traverse_postorder_recursive(root, cbtrav, &index, outbuffer);
  return index;
}

static int rnode_height_recursive(rnode_t * node)
{
  if (!node) return 1;

  int a = rnode_height_recursive(node->left);
  int b = rnode_height_recursive(node->right);

  return MAX(a,b)+1;
}


int rnode_height(rnode_t * root)
{
  return rnode_height_recursive(root);
}

static void rnode_query_tipnodes_recursive(rnode_t * node,
                                           rnode_t ** node_list,
                                           int * index)
{
  if (!node) return;

  if (!node->left)
  {
    node_list[*index] = node;
    *index = *index + 1;
    return;
  }

  rnode_query_tipnodes_recursive(node->left,  node_list, index);
  rnode_query_tipnodes_recursive(node->right, node_list, index);
}

int rnode_query_tipnodes(rnode_t * root,
                         rnode_t ** node_list)
{
  int index = 0;

  if (!root) return 0;
  if (!root->left)
  {
    node_list[index++] = root;
    return index;
  }

  rnode_query_tipnodes_recursive(root->left,  node_list, &index);
  rnode_query_tipnodes_recursive(root->right, node_list, &index);

  return index;
}

static void rnode_query_innernodes_recursive(rnode_t * root,
                                             rnode_t ** node_list,
                                             int * index)
{
  if (!root) return;
  if (!root->left) return;

  /* postorder traversal */

  rnode_query_innernodes_recursive(root->left,  node_list, index);
  rnode_query_innernodes_recursive(root->right, node_list, index);

  node_list[*index] = root;
  *index = *index + 1;
  return;
}

int rnode_query_innernodes(rnode_t * root,
                           rnode_t ** node_list)
{
  int index = 0;

  if (!root) return 0;
  if (!root->left) return 0;

  rnode_query_innernodes_recursive(root->left,  node_list, &index);
  rnode_query_innernodes_recursive(root->right, node_list, &index);

  node_list[index++] = root;

  return index;
}

void rnode_reset_info(rnode_t * root)
{
  if (!root->left)
  {
    root->leaves = 1;
    root->edge_count = 0;
    root->edgelen_sum = 0;
    root->max_species_count = 1;
    return;
  }

  rnode_reset_info(root->left);
  rnode_reset_info(root->right);

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

  root->max_species_count = 1;
  if (root->edge_count > 0)
    root->max_species_count = root->left->max_species_count +
                              root->right->max_species_count;
}

void rnode_print_tips(rnode_t * node, FILE * out)
{
  if (node->left)  rnode_print_tips(node->left,out);
  if (node->right) rnode_print_tips(node->right,out);

  if (!node->left && !node->right)
    fprintf(out, "%s\n", node->label);
}


rnode_t * rnode_clone(rnode_t * node, rnode_t * parent)
{
  if (!node) return NULL;

  /* clone node */
  rnode_t * clone = (rnode_t *)xcalloc(1,sizeof(rnode_t));
  memcpy(clone,node,sizeof(rnode_t));
  clone->parent = parent;
  clone->data = NULL;

  if (node->label)
    clone->label = xstrdup(node->label);

  /* clone the two subtrees */
  clone->left  = rnode_clone(node->left, clone);
  clone->right = rnode_clone(node->right, clone);

  return clone;
}

static rnode_t ** rnode_tipstring_nodes(rnode_t * root,
                                        char * tipstring,
                                        unsigned int * tiplist_count)
{
  size_t i;
  unsigned int k;
  unsigned int commas_count = 0;

  char * taxon;
  unsigned long taxon_len;

  for (i = 0; i < strlen(tipstring); ++i)
    if (tipstring[i] == ',')
      commas_count++;

  rnode_t ** node_list = (rnode_t **)xmalloc((size_t)(root->leaves) *
                                             sizeof(rnode_t *));
  rnode_query_tipnodes(root, node_list);

  rnode_t ** out_node_list = (rnode_t **)xmalloc((size_t)(commas_count+1) *
                                                 sizeof(rnode_t *));

  /* create a hashtable of tip labels */
  hashtable_t * ht = hashtable_create((unsigned long)(root->leaves));

  for (i = 0; i < (unsigned int)(root->leaves); ++i)
  {
    pair_t * pair = (pair_t *)xmalloc(sizeof(pair_t));
    pair->label = node_list[i]->label;
    pair->index = i;

    if (!hashtable_insert(ht,
                          (void *)pair,
                          hash_fnv(node_list[i]->label),
                          hashtable_paircmp))
      fatal("Duplicate taxon (%s)\n", node_list[i]->label);
  }

  char * s = tipstring;

  k = 0;
  while (*s)
  {
    /* get next tip */
    taxon_len = strcspn(s, ",");
    if (!taxon_len)
      fatal("Erroneous prune list format (double comma)/taxon missing");

    taxon = xstrndup(s, taxon_len);

    /* search tip in hash table */
    pair_t * query = hashtable_find(ht,
                                    taxon,
                                    hash_fnv(taxon),
                                    hashtable_paircmp);

    if (!query)
      fatal("Taxon %s in does not appear in the tree", taxon);

    /* store pointer in output list */
    out_node_list[k++] = node_list[query->index];

    /* free tip label, and move to the beginning of next tip if available */
    free(taxon);
    s += taxon_len;
    if (*s == ',')
      s += 1;
  }

  /* kill the hash table */
  hashtable_destroy(ht,free);

  free(node_list);

  /* return number of tips in the list */
  *tiplist_count = commas_count + 1;

  /* return tip node list */
  return out_node_list;
}

/* fill path with nodes of the path tip to root */
static void fill_path(rnode_t ** path, int * path_len, rnode_t * tip)
{
  int i = 0;

  while (tip)
  {
    path[i++] = tip;
    tip = tip->parent;
  }

  *path_len = i;
}

rnode_t * rnode_lca(rnode_t * root,
                    rnode_t ** tip_nodes,
                    unsigned int count)
{
  unsigned int i;
  rnode_t *** path;

  assert(count >= 2);

  /* allocate path arrays for count tip nodes */
  path = (rnode_t ***)xmalloc((size_t)count *
                                  sizeof(rnode_t **));
  int * path_len = (int *)xmalloc((size_t)count * sizeof(int));

  /* for each tip node fill corresponding path array with all nodes
     in the path to the root node and store the length of the path  */
  for (i = 0; i < count; ++i)
  {
    path[i] = (rnode_t **)xmalloc((size_t)(rnode_height(root)) *
                                  sizeof(rnode_t *));

    fill_path(path[i], &(path_len[i]), tip_nodes[i]);
  }

  /* find the LCA using a breadth-first-search traversal starting from the root.
     Since all paths start at the root, the LCA is the parent of nodes that
     differ in the paths when encountered for the first time */
  rnode_t * lca = NULL;
  while (!lca)
  {
    for (i = 0; i < count; ++i)
      --path_len[i];

    for (i = 1; i < count; ++i)
    {
      if (path[i-1][path_len[i-1]] != path[i][path_len[i]])
      {
        lca = path[i][path_len[i]+1];
        break;
      }
    }
  }

  /* free allocated memory */
  for (i = 0; i < count; ++i)
    free(path[i]);
  free(path);
  free(path_len);

  return lca;
}

rnode_t * rtree_outgroup_lca(rtree_t * tree)
{
  rnode_t * root = tree->root;
  unsigned int og_tips_count;
  rnode_t * og_root;
  rnode_t ** og_tips;


  og_tips = rnode_tipstring_nodes(root,
                                  opt_outgroup,
                                  &og_tips_count);

  if (og_tips_count > 1)
    og_root = rnode_lca(root, og_tips, og_tips_count);
  else og_root = og_tips[0];

  free(og_tips);

  return og_root;
}

int rtree_crop(rtree_t * tree, rnode_t * crop_root)
{
  rnode_t * root = tree->root;

  /* check if the selected subtree can be cropped */
  if (root->leaves - crop_root->leaves < 2)
    return -1;

  rnode_t * new_root;

  /* subtree can be cropped, distinguish between two cases: */

  if (crop_root->parent == root)
  {

    /* Case 1:

          root
         *
        / \                               A
     A *   * crop_root     ---->          *
          / \
         *   *

       in this case the subtree rooted at crop_root is cropped, the root node is
       eliminated and subtree rooted at A becomes the new tree
    */

    if (root->left == crop_root)
    {
      new_root = root->right;
      root->right = NULL;
    }
    else
    {
      new_root = root->left;
      root->left = NULL;
    }

    rnode_destroy(root);

    new_root->parent = NULL;
    rnode_reset_info(new_root);
  }
  else
  {

    /* Case 2:

          root
         *
        / \
     A *   -
            \                               root
             * B             ---->         *
            / \                           / \
         C *   * crop_root             A *   -
              / \                             \
             *   *                             * C

       in this case the subtree rooted at crop_root is cropped, the root node is
       eliminated and subtree rooted at A becomes the new tree
    */

    rnode_t * b = crop_root->parent;
    rnode_t * c;

    /* get C and break the link between B and C */
    if (b->left == crop_root)
    {
      c = b->right;
      b->right = NULL;
    }
    else
    {
      c = b->left;
      b->left = NULL;
    }

    /* link the parent of B with C from both directions */
    c->parent = b->parent;
    if (b->parent->left == b)
      b->parent->left = c;
    else
      b->parent->right = c;

    c->length += b->length;

    rnode_destroy(b);
    rnode_reset_info(root);

    new_root = root;
  }

  /* update the tree container */
  free(tree->nodes);
  tree->root = new_root;
  tree->tip_count = (unsigned int)new_root->leaves;
  tree->inner_count = tree->tip_count - 1;
  tree->edge_count = (unsigned int)new_root->edge_count;

  tree->nodes = (rnode_t **)xmalloc(
      (size_t)(tree->tip_count + tree->inner_count) * sizeof(rnode_t *));
  rnode_query_tipnodes(new_root, tree->nodes);
  rnode_query_innernodes(new_root, tree->nodes + tree->tip_count);

  return 0;
}
