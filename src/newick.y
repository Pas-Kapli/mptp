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

extern int yylex();
extern FILE * yyin;
extern void yylex_destroy();

void yyerror(rtree_t * rtree, const char * s) 
{
  fprintf(stderr, "%s.\n", s);
}

%}


%union
{
  char * s;
  char * d;
  struct rtree_noderec * rtree;
}

%error-verbose
%parse-param {struct rtree_noderec * rtree}
%destructor { yy_dealloc_rtree($$); } subtree

%token OPAR
%token CPAR
%token COMMA
%token COLON SEMICOLON 
%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label
%type<d> number optional_length
%type<rtree> subtree
%start input
%%

input: OPAR subtree COMMA subtree CPAR optional_label optional_length SEMICOLON
{
  rtree->left   = $2;
  rtree->right  = $4;
  rtree->label  = $6;
  rtree->length = $7 ? atof($7) : 0;
  rtree->leaves = $2->leaves + $4->leaves;
  rtree->parent = NULL;
  free($7);

  rtree->left->parent  = rtree; 
  rtree->right->parent = rtree; 
};

subtree: OPAR subtree COMMA subtree CPAR optional_label optional_length
{
  $$ = yy_create_rtree();
  $$->left   = $2;
  $$->right  = $4;
  $$->label  = $6;
  $$->length = $7 ? atof($7) : 0;
  $$->leaves = $2->leaves + $4->leaves;
  free($7);

  $$->left->parent  = $$;
  $$->right->parent = $$;
}
       | label optional_length
{
  $$ = yy_create_rtree();
  $$->label  = $1;
  $$->length = $2 ? atof($2) : 0;
  $$->left   = NULL;
  $$->right  = NULL;
  $$->leaves = 1;
  free($2);
};

 
optional_label:  { $$ = NULL;} | label  {$$ = $1;};
optional_length: { $$ = NULL;} | COLON number {$$ = $2;};
label: STRING    { $$=$1;};
number: NUMBER   { $$=$1;};

%%

rtree_t * yy_parse_rtree(const char * filename)
{
  struct rtree_noderec * rtree;

  rtree = yy_create_rtree();

  yyin = fopen(filename, "r");
  if (!yyin)
  {
    fatal("Cannot open file %s", filename);
  }
  else if (yyparse(rtree))
  {
    fatal("Cannot parse tree file %s (maybe non-binary?)", filename);
  }
  
  if (yyin) fclose(yyin);

  yylex_destroy();

  return rtree;
}
