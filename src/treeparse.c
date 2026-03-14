/*
    Copyright (C) 2015-2026 Tomas Flouri

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

/* -- Token types ---------------------------------------------------------- */

#define TOKEN_OPAR      0
#define TOKEN_CPAR      1
#define TOKEN_COMMA     2
#define TOKEN_COLON     3
#define TOKEN_SEMICOLON 4
#define TOKEN_STRING    5

typedef struct ltoken_s
{
  int type;
  char * data;      /* for TOKEN_STRING only */
} ltoken_t;

/* -- Character classification --------------------------------------------- */

/* Characters that cannot start or appear in an unquoted label.
   In newick: whitespace, structural chars ( ) , : ; [ ] and quotes.
   We match the flex lexer's character classes exactly. */
static int is_label_first(int c)
{
  /* the flex pattern is:
     [^ \'\",\(\):;\[\]\t\n\r]  for first char
     So these chars CANNOT start a label. */
  if (c <= 0) return 0;
  switch (c)
  {
    case ' ': case '\t': case '\n': case '\r':
    case '\'': case '"': case ',':
    case '(': case ')': case ':': case ';':
    case '[': case ']':
      return 0;
    default:
      return 1;
  }
}

static int is_label_rest(int c)
{
  /* the flex pattern continuation is:
     [^ \t\n\r\)\(\[\]\,:;]*
     So these chars break a label. */
  if (c <= 0) return 0;
  switch (c)
  {
    case ' ': case '\t': case '\n': case '\r':
    case '(': case ')': case '[': case ']':
    case ',': case ':': case ';':
      return 0;
    default:
      return 1;
  }
}

/* -- Token deallocation --------------------------------------------------- */

static void token_dealloc(void * data)
{
  ltoken_t * token = (ltoken_t *)data;
  if (token->data)
    free(token->data);
  free(token);
}

/* -- Tokenizer ------------------------------------------------------------ */

static list_t * tokenize(const char * s)
{
  list_t * tokens = (list_t *)xcalloc(1, sizeof(list_t));

  while (*s)
  {
    /* skip whitespace */
    if (*s == ' ' || *s == '\t' || *s == '\n' || *s == '\r')
    {
      s++;
      continue;
    }

    ltoken_t * token = (ltoken_t *)xcalloc(1, sizeof(ltoken_t));

    switch (*s)
    {
      case '(':
        token->type = TOKEN_OPAR;
        s++;
        break;
      case ')':
        token->type = TOKEN_CPAR;
        s++;
        break;
      case ',':
        token->type = TOKEN_COMMA;
        s++;
        break;
      case ':':
        token->type = TOKEN_COLON;
        s++;
        break;
      case ';':
        token->type = TOKEN_SEMICOLON;
        s++;
        break;

      case '\'':
      case '"':
      {
        /* quoted string — handle escape sequences matching the flex lexer */
        char quote = *s;
        s++;

        /* accumulate into a dynamically grown buffer */
        size_t cap = 64;
        size_t len = 0;
        char * buf = (char *)xmalloc(cap);

        while (*s && *s != quote)
        {
          if (*s == '\\')
          {
            s++;
            if (!*s)
              fatal("Unterminated escape sequence in quoted label");

            /* match escape behavior from lex_rtree.l */
            if (*s == quote)
            {
              /* escaped quote: store the escape + quote like the flex lexer */
              if (len + 2 >= cap) { cap *= 2; buf = xrealloc(buf, cap); }
              buf[len++] = '\\';
              buf[len++] = *s;
              s++;
            }
            else if (*s == 'n')
            {
              if (len + 2 >= cap) { cap *= 2; buf = xrealloc(buf, cap); }
              buf[len++] = '\\';
              buf[len++] = 'n';
              s++;
            }
            else if (*s == 't')
            {
              if (len + 2 >= cap) { cap *= 2; buf = xrealloc(buf, cap); }
              buf[len++] = '\\';
              buf[len++] = 't';
              s++;
            }
            else if (*s == '\\')
            {
              if (len + 2 >= cap) { cap *= 2; buf = xrealloc(buf, cap); }
              buf[len++] = '\\';
              buf[len++] = '\\';
              s++;
            }
            else
            {
              /* bare backslash (not followed by a recognized escape) */
              if (len + 1 >= cap) { cap *= 2; buf = xrealloc(buf, cap); }
              buf[len++] = '\\';
              /* don't consume *s — it's the next char */
            }
          }
          else
          {
            /* in a single-quoted region, an unescaped double-quote
               is stored literally, and vice versa — matching flex */
            if (len + 1 >= cap) { cap *= 2; buf = xrealloc(buf, cap); }
            buf[len++] = *s;
            s++;
          }
        }

        if (*s != quote)
          fatal("Unterminated quoted label in newick string");
        s++; /* consume closing quote */

        buf[len] = '\0';
        token->type = TOKEN_STRING;
        token->data = buf;
        break;
      }

      case '[':
        /* skip newick comments [...]  */
        s++;
        while (*s && *s != ']')
          s++;
        if (*s == ']')
          s++;
        free(token);
        continue;

      default:
        /* unquoted label or number */
        if (is_label_first((unsigned char)*s))
        {
          const char * start = s;
          s++;
          while (*s && is_label_rest((unsigned char)*s))
            s++;
          token->type = TOKEN_STRING;
          token->data = xstrndup(start, (size_t)(s - start));
        }
        else
        {
          fatal("Syntax error: unexpected character '%c' in newick string", *s);
        }
        break;
    }

    list_append(tokens, token);
  }

  return tokens;
}

/* -- Recursive descent parser --------------------------------------------- */

/*
   Grammar:
     tree     -> subtree ';'
     subtree  -> '(' children ')' label? length?   (inner node)
              |  label length?                      (leaf)
     children -> subtree (',' subtree)*
     label    -> STRING
     length   -> ':' STRING
*/

typedef struct parser_s
{
  list_item_t * current;
} parser_t;

static ltoken_t * peek(parser_t * p)
{
  if (!p->current)
    return NULL;
  return (ltoken_t *)(p->current->data);
}

static ltoken_t * advance(parser_t * p)
{
  if (!p->current)
    return NULL;
  ltoken_t * tok = (ltoken_t *)(p->current->data);
  p->current = p->current->next;
  return tok;
}

static int check(parser_t * p, int type)
{
  ltoken_t * tok = peek(p);
  return (tok && tok->type == type);
}

static void expect(parser_t * p, int type)
{
  ltoken_t * tok = advance(p);
  if (!tok || tok->type != type)
  {
    const char * names[] = {"(", ")", ",", ":", ";", "STRING"};
    if (tok)
      fatal("Newick parse error: expected '%s', got '%s'",
            names[type], names[tok->type]);
    else
      fatal("Newick parse error: expected '%s', got end of input",
            names[type]);
  }
}

static node_t * parse_subtree(parser_t * p);

static node_t * parse_subtree(parser_t * p)
{
  node_t * node = (node_t *)xcalloc(1, sizeof(node_t));

  if (check(p, TOKEN_OPAR))
  {
    /* inner node: '(' children ')' optional_label optional_length */
    advance(p); /* consume '(' */

    /* parse first child */
    int alloc = 4;
    node->children = (node_t **)xmalloc((size_t)alloc * sizeof(node_t *));
    node->children[0] = parse_subtree(p);
    node->children[0]->parent = node;
    node->children_count = 1;

    /* parse remaining children separated by ',' */
    while (check(p, TOKEN_COMMA))
    {
      advance(p); /* consume ',' */

      if (node->children_count >= alloc)
      {
        alloc *= 2;
        node->children = (node_t **)xrealloc(node->children,
                                              (size_t)alloc * sizeof(node_t *));
      }
      node->children[node->children_count] = parse_subtree(p);
      node->children[node->children_count]->parent = node;
      node->children_count++;
    }

    expect(p, TOKEN_CPAR);

    /* optional label */
    if (check(p, TOKEN_STRING))
    {
      ltoken_t * tok = advance(p);
      node->label = tok->data;
      tok->data = NULL; /* transfer ownership */
    }

    /* optional length */
    if (check(p, TOKEN_COLON))
    {
      advance(p); /* consume ':' */
      if (!check(p, TOKEN_STRING))
        fatal("Newick parse error: expected branch length after ':'");
      ltoken_t * tok = advance(p);
      node->length = atof(tok->data);
    }
  }
  else if (check(p, TOKEN_STRING))
  {
    /* leaf: label optional_length */
    ltoken_t * tok = advance(p);
    node->label = tok->data;
    tok->data = NULL; /* transfer ownership */

    /* optional length */
    if (check(p, TOKEN_COLON))
    {
      advance(p);
      if (!check(p, TOKEN_STRING))
        fatal("Newick parse error: expected branch length after ':'");
      tok = advance(p);
      node->length = atof(tok->data);
    }

    node->leaves = 1;
  }
  else
  {
    /* empty label internal node or empty (degenerate but handle gracefully) */
    if (check(p, TOKEN_COLON))
    {
      advance(p);
      if (!check(p, TOKEN_STRING))
        fatal("Newick parse error: expected branch length after ':'");
      ltoken_t * tok = advance(p);
      node->length = atof(tok->data);
    }
  }

  return node;
}

/* -- Count leaves in ntree ------------------------------------------------ */

static int count_leaves(node_t * node)
{
  if (node->children_count == 0)
    return 1;

  int sum = 0;
  int i;
  for (i = 0; i < node->children_count; i++)
    sum += count_leaves(node->children[i]);
  return sum;
}

/* -- Build ntree from parse ----------------------------------------------- */

static void collect_nodes(node_t * node, node_t ** leaves, int * leaf_idx,
                          node_t ** inner, int * inner_idx)
{
  int i;
  if (node->children_count == 0)
  {
    leaves[*leaf_idx] = node;
    (*leaf_idx)++;
    return;
  }

  for (i = 0; i < node->children_count; i++)
    collect_nodes(node->children[i], leaves, leaf_idx, inner, inner_idx);

  inner[*inner_idx] = node;
  (*inner_idx)++;
}

/* -- Public API: parse newick string into ntree_t ------------------------- */

ntree_t * parse_newick_string(const char * newick)
{
  list_t * tokens = tokenize(newick);

  if (tokens->count == 0)
  {
    list_clear(tokens, token_dealloc);
    free(tokens);
    return NULL;
  }

  parser_t parser;
  parser.current = tokens->head;

  node_t * root = parse_subtree(&parser);

  /* expect semicolon */
  if (check(&parser, TOKEN_SEMICOLON))
    advance(&parser);

  /* compute leaf counts */
  root->leaves = count_leaves(root);

  /* set leaves recursively */
  /* We'll use a simple recursive pass */
  /* (leaves already set for leaf nodes during parse, inner needs computation) */

  /* build ntree_t */
  ntree_t * tree = (ntree_t *)xcalloc(1, sizeof(ntree_t));
  tree->root = root;
  tree->tip_count = root->leaves;

  /* count inner nodes */
  int inner_count = 0;
  {
    /* simple recursive count */
    /* inner nodes = total nodes - tip count */
    /* total_nodes = tip_count + inner_count  */
    /* but we need to count actual inner nodes */
    node_t ** stack = (node_t **)xmalloc((size_t)(2 * tree->tip_count) * sizeof(node_t *));
    int sp = 0;
    stack[sp++] = root;
    int total = 0;
    while (sp > 0)
    {
      node_t * n = stack[--sp];
      total++;
      int i;
      for (i = 0; i < n->children_count; i++)
        stack[sp++] = n->children[i];
    }
    free(stack);
    inner_count = total - tree->tip_count;
  }
  tree->inner_count = inner_count;

  /* populate node arrays */
  tree->leaves_list = (node_t **)xmalloc((size_t)tree->tip_count * sizeof(node_t *));
  tree->inner_list = (node_t **)xmalloc((size_t)tree->inner_count * sizeof(node_t *));

  int leaf_idx = 0;
  int inner_idx = 0;
  collect_nodes(root, tree->leaves_list, &leaf_idx,
                tree->inner_list, &inner_idx);

  /* clean up token list (data was transferred, so just free token structs) */
  list_clear(tokens, token_dealloc);
  free(tokens);

  return tree;
}

/* -- Classify and validate tree ------------------------------------------- */

typedef enum
{
  TREE_ROOTED_BINARY,
  TREE_UNROOTED_BINARY,
  TREE_MULTIFURCATING,
  TREE_INVALID
} tree_class_t;

static tree_class_t classify_tree(ntree_t * tree)
{
  int i;
  int root_children = tree->root->children_count;

  if (root_children < 2)
    return TREE_INVALID;

  /* check all inner nodes (except root) for binary structure */
  for (i = 0; i < tree->inner_count; i++)
  {
    node_t * n = tree->inner_list[i];
    if (n == tree->root)
      continue;
    if (n->children_count != 2)
      return TREE_MULTIFURCATING;
  }

  if (root_children == 2)
    return TREE_ROOTED_BINARY;
  else if (root_children == 3)
  {
    /* verify root also has binary children (already checked above) */
    return TREE_UNROOTED_BINARY;
  }
  else
    return TREE_MULTIFURCATING;
}

/* -- Convert ntree to rnode (rooted binary) ------------------------------- */

static rnode_t * ntree_to_rnode(node_t * nnode)
{
  rnode_t * rnode = (rnode_t *)xcalloc(1, sizeof(rnode_t));

  rnode->event = EVENT_COALESCENT;
  rnode->mark = 0;
  rnode->data = NULL;

  if (nnode->label)
    rnode->label = xstrdup(nnode->label);
  else
    rnode->label = NULL;
  rnode->length = nnode->length;

  if (nnode->children_count == 0)
  {
    /* leaf */
    rnode->left = NULL;
    rnode->right = NULL;
    rnode->leaves = 1;
    rnode->edge_count = 0;
    rnode->edgelen_sum = 0;
    rnode->max_species_count = 1;
    return rnode;
  }

  /* inner node with exactly 2 children */
  rnode->left = ntree_to_rnode(nnode->children[0]);
  rnode->right = ntree_to_rnode(nnode->children[1]);

  rnode->left->parent = rnode;
  rnode->right->parent = rnode;

  rnode->leaves = rnode->left->leaves + rnode->right->leaves;

  rnode->edge_count = rnode->left->edge_count + rnode->right->edge_count;
  rnode->edgelen_sum = rnode->left->edgelen_sum + rnode->right->edgelen_sum;

  if (rnode->left->length > opt_minbr)
  {
    rnode->edge_count++;
    rnode->edgelen_sum += rnode->left->length;
  }
  if (rnode->right->length > opt_minbr)
  {
    rnode->edge_count++;
    rnode->edgelen_sum += rnode->right->length;
  }

  rnode->max_species_count = 1;
  if (rnode->edge_count > 0)
    rnode->max_species_count = rnode->left->max_species_count +
                               rnode->right->max_species_count;

  return rnode;
}

static rtree_t * make_rtree(rnode_t * root)
{
  int tip_count = root->leaves;
  int inner_count = tip_count - 1;

  rtree_t * tree = (rtree_t *)xcalloc(1, sizeof(rtree_t));
  tree->root = root;
  tree->tip_count = (unsigned int)tip_count;
  tree->inner_count = (unsigned int)inner_count;
  tree->edge_count = (unsigned int)root->edge_count;

  /* populate nodes array: tips [0..tip_count-1], inner [tip_count..] */
  tree->nodes = (rnode_t **)xmalloc((size_t)(tip_count + inner_count) *
                                    sizeof(rnode_t *));

  rnode_t ** tip_list = tree->nodes;
  rnode_t ** inner_list = tree->nodes + tip_count;

  int tip_idx = 0;
  int inner_idx = 0;

  rnode_query_tipnodes(root, tip_list);
  rnode_query_innernodes(root, inner_list);

  (void)tip_idx;
  (void)inner_idx;

  return tree;
}

/* -- Unrooted to rooted conversion ---------------------------------------- */

/* Find outgroup among ntree children (for unrooted tree) */

static node_t * ntree_find_longest_branchtip(ntree_t * tree)
{
  int i;
  double best_length = -1;
  node_t * best = NULL;

  for (i = 0; i < tree->tip_count; i++)
  {
    if (tree->leaves_list[i]->length > best_length)
    {
      best_length = tree->leaves_list[i]->length;
      best = tree->leaves_list[i];
    }
  }

  return best;
}

/* Find LCA of outgroup tips in the ntree using hash table (matching utree.c logic) */

static node_t * ntree_find_outgroup_lca(ntree_t * tree)
{
  unsigned int i;
  unsigned int k;
  unsigned int commas_count = 0;
  char * taxon;
  size_t taxon_len;

  /* count commas to determine number of outgroup taxa */
  for (i = 0; i < strlen(opt_outgroup); i++)
    if (opt_outgroup[i] == ',')
      commas_count++;

  unsigned int og_count = commas_count + 1;

  /* create hash table of all tip labels */
  hashtable_t * ht = hashtable_create((unsigned long)tree->tip_count);
  for (i = 0; i < (unsigned int)tree->tip_count; i++)
  {
    pair_t * pair = (pair_t *)xmalloc(sizeof(pair_t));
    pair->label = tree->leaves_list[i]->label;
    pair->index = i;

    if (!hashtable_insert(ht, (void *)pair,
                          hash_fnv(tree->leaves_list[i]->label),
                          hashtable_paircmp))
      fatal("Duplicate taxon (%s)\n", tree->leaves_list[i]->label);
  }

  /* find outgroup tips */
  node_t ** og_tips = (node_t **)xmalloc((size_t)og_count * sizeof(node_t *));
  char * s = opt_outgroup;
  k = 0;
  while (*s)
  {
    taxon_len = strcspn(s, ",");
    if (!taxon_len)
      fatal("Erroneous outgroup format (double comma)/taxon missing");

    taxon = xstrndup(s, taxon_len);

    pair_t * query = hashtable_find(ht, taxon, hash_fnv(taxon),
                                    hashtable_paircmp);
    if (!query)
      fatal("Taxon %s does not appear in the tree", taxon);

    og_tips[k++] = tree->leaves_list[query->index];

    free(taxon);
    s += taxon_len;
    if (*s == ',')
      s += 1;
  }

  hashtable_destroy(ht, free);

  /* if single outgroup tip, return it directly */
  if (og_count == 1)
  {
    node_t * result = og_tips[0];
    free(og_tips);
    return result;
  }

  /* For multiple outgroup taxa, find their LCA.
     Mark all outgroup tips, then walk up from each to find common ancestor */

  /* Mark outgroup tips */
  for (i = 0; i < og_count; i++)
    og_tips[i]->data = (void *)1;

  /* Walk from first outgroup tip to root, storing path */
  int max_depth = tree->tip_count + tree->inner_count;
  node_t ** path = (node_t **)xmalloc((size_t)max_depth * sizeof(node_t *));
  int path_len = 0;
  node_t * cur = og_tips[0];
  while (cur)
  {
    path[path_len++] = cur;
    cur = cur->parent;
  }

  /* For each path node (from tip to root), check if it's ancestor of all
     outgroup tips. The first one from the bottom (closest to tips) that is
     ancestor of all is the LCA. */

  /* Simple approach: mark path nodes, then for each og tip, walk up and
     find the first marked node */

  /* Better approach: for each candidate from tip toward root, count
     how many og tips are in its subtree */

  /* Simplest correct approach: iterative LCA for the set */
  node_t * lca = og_tips[0];
  for (i = 1; i < og_count; i++)
  {
    /* Find LCA of lca and og_tips[i] */
    /* Mark path from lca to root */
    cur = lca;
    while (cur)
    {
      cur->data = (void *)2;
      cur = cur->parent;
    }
    /* Walk from og_tips[i] to root, find first marked */
    cur = og_tips[i];
    while (cur)
    {
      if (cur->data == (void *)2)
      {
        lca = cur;
        break;
      }
      cur = cur->parent;
    }
    /* Unmark path */
    cur = lca;
    while (cur)
    {
      if (cur->data == (void *)2)
        cur->data = NULL;
      cur = cur->parent;
    }
  }

  /* Clear all marks */
  for (i = 0; i < og_count; i++)
    og_tips[i]->data = NULL;

  free(og_tips);
  free(path);

  /* Validate: the LCA must have exactly the outgroup tips in its subtree.
     Count leaves under lca and compare with og_count */
  int lca_leaves = count_leaves(lca);
  if ((unsigned int)lca_leaves != og_count)
  {
    return NULL; /* not a valid monophyletic outgroup */
  }

  return lca;
}

/* Convert an unrooted ntree (root with 3 children) to rooted rtree.
   Root at the outgroup edge: the outgroup becomes one side, everything
   else becomes the other side, with branch length split in half. */

static rtree_t * ntree_unrooted_to_rtree(ntree_t * tree)
{
  node_t * og_node;

  if (!opt_outgroup)
  {
    og_node = ntree_find_longest_branchtip(tree);
    if (!og_node)
      fatal("Failed to find outgroup tip");
    fprintf(stdout,
            "Selected %s as outgroup based on longest tip-branch criterion\n",
            og_node->label);
  }
  else
  {
    og_node = ntree_find_outgroup_lca(tree);
    if (!og_node)
      fatal("Outgroup must be a single tip or a list of all tips of a subtree");
  }

  /* Handle crop */
  if (opt_crop)
  {
    /* We need to build a rooted tree from the non-outgroup side.
       Find the parent of the outgroup node. The other children of
       that parent become the rooted tree. */

    node_t * parent = og_node->parent;
    if (!parent)
      fatal("Cannot crop: outgroup is the root");

    /* For unrooted trees, the outgroup should be a child of the root
       (which has 3 children). If it's deeper, we need to handle that. */

    /* If outgroup is a direct child of root (the common case for unrooted):
       root has 3 children, remove the outgroup one, root the remaining 2 */
    if (parent == tree->root)
    {
      /* find the two non-outgroup children */
      node_t * others[2];
      int oc = 0;
      int i;
      for (i = 0; i < parent->children_count; i++)
      {
        if (parent->children[i] != og_node)
        {
          if (oc >= 2)
            fatal("Unexpected number of children when cropping");
          others[oc++] = parent->children[i];
        }
      }
      if (oc != 2)
        fatal("Expected exactly 2 remaining children after crop");

      /* Build rooted tree from these two subtrees */
      rnode_t * root = (rnode_t *)xcalloc(1, sizeof(rnode_t));
      root->left = ntree_to_rnode(others[0]);
      root->right = ntree_to_rnode(others[1]);
      root->left->parent = root;
      root->right->parent = root;
      root->parent = NULL;
      root->length = 0;
      root->label = NULL;
      root->data = NULL;
      root->mark = 0;
      root->event = EVENT_COALESCENT;

      rnode_reset_info(root);
      return make_rtree(root);
    }
    else
    {
      /* Outgroup is deeper in the tree. Remove the outgroup subtree by
         bypassing its parent: replace parent with sibling in grandparent's
         children list, then build the rooted tree directly from the modified
         root's 3 children (which now excludes the outgroup).

         We must NOT use goto convert_unrooted here because that path uses
         og_node for rerooting, but og_node's parent linkage is now invalid
         after the bypass surgery. Instead, build the rooted tree directly
         from the remaining 3-child root. */

      /* parent has children_count == 2 (binary inner node) */
      if (parent->children_count != 2)
        fatal("Cannot crop: outgroup parent is not binary");

      node_t * sibling = NULL;
      int i;
      for (i = 0; i < parent->children_count; i++)
        if (parent->children[i] != og_node)
          sibling = parent->children[i];

      if (!sibling)
        fatal("Cannot find sibling of outgroup");

      /* Replace parent in grandparent's children with sibling */
      /* and add parent's branch length to sibling's */
      sibling->length += parent->length;
      sibling->parent = parent->parent;

      if (parent->parent)
      {
        for (i = 0; i < parent->parent->children_count; i++)
        {
          if (parent->parent->children[i] == parent)
          {
            parent->parent->children[i] = sibling;
            break;
          }
        }
      }

      /* The root still has 3 children but one subtree has been modified
         (the parent node was bypassed). Collect the root's 3 children
         and build a rooted tree from them directly. */

      node_t * others[3];
      int oc = 0;
      for (i = 0; i < tree->root->children_count; i++)
      {
        if (oc >= 3)
          fatal("Root has more than 3 children after crop");
        others[oc++] = tree->root->children[i];
      }

      if (oc != 3)
        fatal("Expected 3 children at unrooted root after crop");

      /* Pick the longest tip branch as outgroup for rooting.
         Use the first child as default, then find the longest. */
      node_t * new_og = others[0];
      double max_len = others[0]->length;
      for (i = 1; i < 3; i++)
      {
        if (others[i]->length > max_len)
        {
          max_len = others[i]->length;
          new_og = others[i];
        }
      }

      /* Build rooted tree: new_og on left, other 2 form right subtree */
      rnode_t * root = (rnode_t *)xcalloc(1, sizeof(rnode_t));
      root->left = ntree_to_rnode(new_og);

      node_t * remaining[2];
      int rc = 0;
      for (i = 0; i < 3; i++)
        if (others[i] != new_og)
          remaining[rc++] = others[i];

      rnode_t * inner = (rnode_t *)xcalloc(1, sizeof(rnode_t));
      inner->left = ntree_to_rnode(remaining[0]);
      inner->right = ntree_to_rnode(remaining[1]);
      inner->left->parent = inner;
      inner->right->parent = inner;
      inner->label = NULL;
      inner->length = 0;
      inner->event = EVENT_COALESCENT;
      inner->data = NULL;
      inner->mark = 0;

      root->right = inner;
      root->left->parent = root;
      root->right->parent = root;

      /* Split the outgroup branch length in half between the two sides */
      root->left->length /= 2;
      root->right->length = root->left->length;

      root->parent = NULL;
      root->length = 0;
      root->label = NULL;
      root->event = EVENT_COALESCENT;
      root->data = NULL;
      root->mark = 0;

      rnode_reset_info(root);
      return make_rtree(root);
    }
  }

  /* Standard unrooted to rooted conversion:
     Find which side of the root the outgroup is on.
     Root at the edge connecting outgroup to the rest. */

  /* If outgroup is a direct child of root */
  if (og_node->parent == tree->root)
  {
    /* Root at the edge from og_node to root.
       og_node becomes one child, the other children form the other subtree */
    rnode_t * root = (rnode_t *)xcalloc(1, sizeof(rnode_t));

    /* Left child is the outgroup */
    root->left = ntree_to_rnode(og_node);

    /* Right child: if root has exactly 3 children, the other 2 form a subtree */
    node_t * others[2];
    int oc = 0;
    int i;
    for (i = 0; i < tree->root->children_count; i++)
    {
      if (tree->root->children[i] != og_node)
      {
        if (oc >= 2)
          fatal("Root has more than 3 children");
        others[oc++] = tree->root->children[i];
      }
    }

    if (oc == 2)
    {
      /* Create an inner rnode for the two non-outgroup children */
      rnode_t * inner = (rnode_t *)xcalloc(1, sizeof(rnode_t));
      inner->left = ntree_to_rnode(others[0]);
      inner->right = ntree_to_rnode(others[1]);
      inner->left->parent = inner;
      inner->right->parent = inner;
      inner->label = NULL;
      inner->length = 0;
      inner->event = EVENT_COALESCENT;
      inner->data = NULL;
      inner->mark = 0;

      root->right = inner;
    }
    else
    {
      /* shouldn't happen for a valid unrooted binary tree */
      fatal("Unexpected number of non-outgroup children");
    }

    root->left->parent = root;
    root->right->parent = root;

    /* Split the outgroup branch length in half between the two sides */
    root->left->length /= 2;
    root->right->length = root->left->length;

    root->parent = NULL;
    root->length = 0;
    root->label = NULL;
    root->event = EVENT_COALESCENT;
    root->data = NULL;
    root->mark = 0;

    rnode_reset_info(root);
    return make_rtree(root);
  }
  else
  {
    /* Outgroup is deeper. We need to re-root at the edge
       from og_node to og_node->parent.

       Walk from og_node up to root, flipping parent-child relationships
       to create a new root at that edge. */

    /* Build outgroup side (just convert og_node subtree) */
    rnode_t * og_rnode = ntree_to_rnode(og_node);

    /* Build the rest of the tree by converting everything except the
       outgroup subtree. This requires "rerooting" the ntree.

       Strategy: we need to create a rooted tree where:
       - One child is the outgroup
       - Other child is the rest of the tree

       We'll do this by creating a new intermediate structure.
       Walk from og_node->parent to the root, reversing directions. */

    /* Collect the path from og_node->parent to root */
    node_t ** rpath = (node_t **)xmalloc(
        (size_t)(tree->tip_count + tree->inner_count) * sizeof(node_t *));
    int rpath_len = 0;
    node_t * cur2 = og_node->parent;
    while (cur2)
    {
      rpath[rpath_len++] = cur2;
      cur2 = cur2->parent;
    }

    /* Build the non-outgroup side by converting nodes along the path.
       At each step, we take the "other" children (not in the path direction)
       and combine them. */

    /* Start from the root (last in rpath) and work back down */

    /* Build from root down to og_node->parent */
    /* The root in ntree has 3 children. Two of them are "not on the path". */
    /* Actually, one child is on the path (or the path passes through one child). */

    /* Let's build it recursively from the root end. */

    /* The new "rest" subtree is constructed by treating the path
       as a spine and attaching the non-path children. */

    /* For the root (3 children): remove the child that's on the path,
       the other 2 children form a subtree. Then walk down the path,
       at each node, one child is the "continuation" (on path), the
       other child becomes a new subtree we attach. */

    /* Special case: if path has length 1 (og_node->parent == root),
       we already handled this above */

    rnode_t * rest = NULL;

    /* Process root node (has 3 children in unrooted tree) */
    node_t * root_nnode = rpath[rpath_len - 1];
    assert(root_nnode == tree->root);

    /* Find which child of root is on the path */
    node_t * path_child = (rpath_len >= 2) ? rpath[rpath_len - 2] : NULL;

    /* Collect non-path children of root */
    node_t ** root_others = (node_t **)xmalloc(
        (size_t)root_nnode->children_count * sizeof(node_t *));
    int root_others_count = 0;
    int i;
    for (i = 0; i < root_nnode->children_count; i++)
    {
      if (root_nnode->children[i] != path_child)
        root_others[root_others_count++] = root_nnode->children[i];
    }

    if (root_others_count == 2)
    {
      /* Create a subtree from these two */
      rest = (rnode_t *)xcalloc(1, sizeof(rnode_t));
      rest->left = ntree_to_rnode(root_others[0]);
      rest->right = ntree_to_rnode(root_others[1]);
      rest->left->parent = rest;
      rest->right->parent = rest;
      rest->label = NULL;
      rest->event = EVENT_COALESCENT;
      rest->data = NULL;
      rest->mark = 0;

      /* The branch length from root to path_child becomes part of the spine */
      rest->length = path_child ? path_child->length : 0;
    }
    else if (root_others_count == 1)
    {
      rest = ntree_to_rnode(root_others[0]);
      rest->length += path_child ? path_child->length : 0;
    }
    else
    {
      fatal("Unexpected number of root children");
    }
    free(root_others);

    /* Now walk down the path from root toward og_node->parent,
       combining the 'rest' with the sibling at each step */
    for (i = rpath_len - 2; i >= 0; i--)
    {
      node_t * current = rpath[i];
      node_t * next_on_path = (i > 0) ? rpath[i - 1] : og_node;

      /* Find the sibling(s) of next_on_path under current */
      int j;
      for (j = 0; j < current->children_count; j++)
      {
        if (current->children[j] != next_on_path)
        {
          /* Convert this sibling */
          rnode_t * sibling_rnode = ntree_to_rnode(current->children[j]);

          /* Combine rest and sibling into new inner node */
          rnode_t * new_inner = (rnode_t *)xcalloc(1, sizeof(rnode_t));
          new_inner->left = sibling_rnode;
          new_inner->right = rest;
          new_inner->left->parent = new_inner;
          new_inner->right->parent = new_inner;
          new_inner->label = current->label ? xstrdup(current->label) : NULL;
          new_inner->length = next_on_path->length;
          new_inner->event = EVENT_COALESCENT;
          new_inner->data = NULL;
          new_inner->mark = 0;

          rest = new_inner;
        }
      }
    }

    free(rpath);

    /* Now create the final root */
    rnode_t * final_root = (rnode_t *)xcalloc(1, sizeof(rnode_t));
    final_root->left = og_rnode;
    final_root->right = rest;
    final_root->left->parent = final_root;
    final_root->right->parent = final_root;

    /* Split branch length at the outgroup edge in half */
    double total_bl = og_rnode->length;
    final_root->left->length = total_bl / 2;
    final_root->right->length = total_bl / 2;

    final_root->parent = NULL;
    final_root->length = 0;
    final_root->label = NULL;
    final_root->event = EVENT_COALESCENT;
    final_root->data = NULL;
    final_root->mark = 0;

    rnode_reset_info(final_root);
    return make_rtree(final_root);
  }
}

/* -- ntree destruction ---------------------------------------------------- */

static void node_destroy(node_t * node)
{
  int i;
  if (!node) return;

  for (i = 0; i < node->children_count; i++)
    node_destroy(node->children[i]);

  free(node->label);
  free(node->children);
  free(node);
}

void ntree_destroy(ntree_t * tree)
{
  if (!tree) return;

  node_destroy(tree->root);
  free(tree->leaves_list);
  free(tree->inner_list);
  free(tree);
}

/* -- Main entry point ----------------------------------------------------- */

rtree_t * parse_newick(const char * filename)
{
  FILE * fp;
  long filesize;
  char * newick;

  fp = fopen(filename, "r");
  if (!fp)
    fatal("Unable to open file (%s)", filename);

  /* get file size */
  if (fseek(fp, 0, SEEK_END) != 0)
    fatal("Unable to seek in file (%s)", filename);
  filesize = ftell(fp);
  if (filesize < 0)
    fatal("Unable to determine size of file (%s)", filename);
  rewind(fp);

  /* read entire file */
  newick = (char *)xmalloc((size_t)(filesize + 1));
  size_t bytes_read = fread(newick, 1, (size_t)filesize, fp);
  newick[bytes_read] = '\0';
  fclose(fp);

  /* parse newick string */
  ntree_t * ntree = parse_newick_string(newick);
  free(newick);

  if (!ntree)
    fatal("Unable to parse tree from file %s", filename);

  /* classify tree */
  tree_class_t tclass = classify_tree(ntree);

  rtree_t * rtree = NULL;

  switch (tclass)
  {
    case TREE_ROOTED_BINARY:
      if (!opt_quiet)
        fprintf(stdout, "Loaded rooted tree...\n");

      /* convert directly */
      {
        rnode_t * root = ntree_to_rnode(ntree->root);
        root->parent = NULL;
        rtree = make_rtree(root);
      }

      /* handle crop for rooted trees */
      if (opt_crop)
      {
        if (!opt_outgroup)
          fatal("--outgroup must be specified when using --outgroup_crop.");

        rnode_t * og_root = rtree_outgroup_lca(rtree);
        if (rtree_crop(rtree, og_root))
          fatal("Cropping the outgroup leads to less than two tips.");
      }
      break;

    case TREE_UNROOTED_BINARY:
      if (!opt_quiet)
      {
        fprintf(stdout, "Loaded unrooted tree...\n");
        fprintf(stdout, "Converting to rooted tree...\n");
      }
      rtree = ntree_unrooted_to_rtree(ntree);
      break;

    case TREE_MULTIFURCATING:
      ntree_destroy(ntree);
      fatal("Tree contains multifurcations (polytomies). "
            "mptp requires fully resolved (binary) trees.");
      break;

    case TREE_INVALID:
      ntree_destroy(ntree);
      fatal("Tree structure is invalid (root has fewer than 2 children).");
      break;
  }

  ntree_destroy(ntree);
  return rtree;
}
