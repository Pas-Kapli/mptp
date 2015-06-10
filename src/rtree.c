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

rtree_t * yy_create_rtree()
{
  rtree_t * t = xrealloc(0, sizeof(rtree_t));
  memset(t, 0, sizeof(rtree_t));
  return t;
}

void yy_dealloc_rtree(rtree_t * tree)
{
  if (!tree) return;

  free(tree->label);
  yy_dealloc_rtree(tree->left);
  yy_dealloc_rtree(tree->right);
  free(tree);
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

  printf (" %s:%f\n", tree->label, tree->length);

  if (active_node_order[indend_level-1] == 2) 
    active_node_order[indend_level-1] = 0;

  active_node_order[indend_level] = 1;
  print_tree_recurse(tree->left, indend_level+1, active_node_order);
  active_node_order[indend_level] = 2;
  print_tree_recurse(tree->right, indend_level+1, active_node_order);

}

int tree_indend_level(rtree_t * tree, int indend)
{
  if (!tree) return indend;

  int a,b;

  a = tree_indend_level(tree->left, indend+1);
  b = tree_indend_level(tree->right, indend+1);

  return MAX(a,b);
}

void show_ascii_rtree(rtree_t * tree)
{
  int indend_max = tree_indend_level(tree,0);
  int * active_node_order = (int *)malloc((indend_max+1) * sizeof(int));
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  printf (" %s:%f\n", tree->label, tree->length);
  print_tree_recurse(tree->left, 1, active_node_order);
  print_tree_recurse(tree->right, 1, active_node_order);
  free(active_node_order);
}

char * export_newick(rtree_t * root)
{
  char * newick;

  if (!root) return NULL;

  if (!(root->left) || !(root->right))
    asprintf(&newick, "%s:%f", root->label, root->length);
  else
  {
    char * subtree1 = export_newick(root->left);
    char * subtree2 = export_newick(root->right);

    asprintf(&newick, "(%s,%s)%s:%f", subtree1,
                                      subtree2,
                                      root->label ? root->label : "",
                                      root->length);
    free(subtree1);
    free(subtree2);
  }

  return newick;
}
