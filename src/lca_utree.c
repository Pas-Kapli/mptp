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

static utree_t ** path1 = NULL;
static utree_t ** path2 = NULL;
static int path1_len = 0;
static int path2_len = 0;

/* fill path with nodes of the path tip to root */
static void fill_path(utree_t ** path, int * path_len, utree_t * tip)
{
  int i = 0;
  
  assert(!tip->next);

  path[i++] = tip;
  tip = tip->back;
  path[i++] = tip;

  while (1)
  {
    if (tip->next->back->height > tip->height)
    {
      tip = tip->next->back;
    }
    else if (tip->next->next->back->height > tip->height)
    {
      tip = tip->next->next->back;
    }
    else
      break;

    path[i++] = tip;
  }

  *path_len = i;
}

void lca_init(utree_t * root)
{
  /* allocate two paths of maximal height */
  path1 = (utree_t **)xmalloc((size_t)(root->height+1) * sizeof(utree_t *));
  path2 = (utree_t **)xmalloc((size_t)(root->height+1) * sizeof(utree_t *));
}

utree_t * lca_compute(utree_t * tip1, utree_t * tip2)
{
  utree_t * outgroup;

  fill_path(path1, &path1_len, tip1);
  fill_path(path2, &path2_len, tip2);

  assert(path1_len && path2_len);

  while(path1[--path1_len] == path2[--path2_len])
  {
    assert(path1_len && path2_len);
  }

  /* get the correct node in the round-about structure such that the subtrees
     of the two outgroups form a common clade */
  if (path1[path1_len]->next == path2[path2_len])
    outgroup = path1[path1_len]->next->next;
  else
    outgroup = path1[path1_len]->next;

  /* return the LCA */
  return outgroup;
}

void lca_destroy()
{
  free(path1);
  free(path2);
}
