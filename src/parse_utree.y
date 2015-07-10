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
%{
#include "delimit.h"

extern int utree_lex();
extern FILE * utree_in;
extern void utree_lex_destroy();

static int tip_cnt = 0;

static void dealloc_tree_recursive(utree_t * node)
{
  if (!node->next)
  {
    free(node->label);
    free(node);
    return;
  }

  dealloc_tree_recursive(node->next->back);
  dealloc_tree_recursive(node->next->next->back);

  free(node->next->next);
  free(node->next);
  free(node->label);
  free(node);
}

void utree_destroy(utree_t * root)
{
  if (!root) return;
  if (!(root->next))
  {
    free(root->label);
    free(root);
    return;
  }

  if (root->next)
    dealloc_tree_recursive(root->next->back);
  if (root->next->next)
    dealloc_tree_recursive(root->next->next->back);
  if (root->back)
    dealloc_tree_recursive(root->back);

  free(root->label);
  free(root->next->next);
  free(root->next);
  free(root);
}

static void utree_error(utree_t * tree, const char * s) 
{
  fprintf(stderr, "%s.\n", s);
}

%}

%union
{
  char * s;
  char * d;
  struct utree_s * tree;
}

%error-verbose
%parse-param {struct utree_s * tree}
%destructor { utree_destroy($$); } subtree

%token OPAR
%token CPAR
%token COMMA
%token COLON SEMICOLON 
%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label
%type<d> number optional_length
%type<tree> subtree
%start input
%%

input: OPAR subtree COMMA subtree COMMA subtree CPAR optional_label optional_length SEMICOLON
{
  tree->next               = (utree_t *)calloc(1, sizeof(utree_t));

  tree->next->next         = (utree_t *)calloc(1, sizeof(utree_t));
  tree->next->next->next   = tree;


  tree->back               = $2;
  tree->next->back         = $4;
  tree->next->next->back   = $6;

  $2->back                 = tree;
  $4->back                 = tree->next;
  $6->back                 = tree->next->next;

  tree->label              = $8;
  tree->next->label        = $8;
  tree->next->next->label  = $8;

  tree->length             = $2->length;
  tree->next->length       = $4->length;
  tree->next->next->length = $6->length;

  free($9);
};

subtree: OPAR subtree COMMA subtree CPAR optional_label optional_length
{
  $$                     = (utree_t *)calloc(1, sizeof(utree_t));

  $$->next               = (utree_t *)calloc(1, sizeof(utree_t));

  $$->next->next         = (utree_t *)calloc(1, sizeof(utree_t));
  $$->next->next->next   = $$;


  $$->next->back         = $2;
  $$->next->next->back   = $4;

  $2->back               = $$->next;
  $4->back               = $$->next->next;

  $$->label              = $6;
  $$->next->label        = $6;
  $$->next->next->label  = $6;

  $$->length = $7 ? atof($7) : 0;
  free($7);

  $$->next->length       = $2->length;
  $$->next->next->length = $4->length;
}
       | label optional_length
{
  $$ = (utree_t *)calloc(1, sizeof(utree_t));

  $$->label  = $1;
  $$->length = $2 ? atof($2) : 0;
  $$->next   = NULL;
  tip_cnt++;
  free($2);
};

 
optional_label:  { $$ = NULL;} | label  {$$ = $1;};
optional_length: { $$ = NULL;} | COLON number {$$ = $2;};
label: STRING    { $$=$1;} | NUMBER {$$=$1;};
number: NUMBER   { $$=$1;};

%%

utree_t * utree_parse_newick(const char * filename, int * tip_count)
{
  struct utree_s * tree;

  /* reset tip count */
  tip_cnt = 0;

  tree = (utree_t *)calloc(1, sizeof(utree_t));

  utree_in = fopen(filename, "r");
  if (!utree_in)
  {
    utree_destroy(tree);
    snprintf(errmsg, 200, "Unable to open file (%s)", filename);
    return NULL;
  }
  else if (utree_parse(tree))
  {
    utree_destroy(tree);
    tree = NULL;
    fclose(utree_in);
    utree_lex_destroy();
    return NULL;
  }
  
  if (utree_in) fclose(utree_in);

  utree_lex_destroy();

  *tip_count = tip_cnt;
  
  return tree;
}
