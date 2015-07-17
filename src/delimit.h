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

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <x86intrin.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdbool.h>

/* constants */

#define PROG_NAME "delimit"
#define PROG_VERSION "v0.0.0"

#ifdef __APPLE__
#define PROG_ARCH "macosx_x86_64"
#else
#define PROG_ARCH "linux_x86_64"
#endif

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;

typedef struct utree_s
{
  char * label;
  double length;
  int height;
  struct utree_s * next;
  struct utree_s * back;

  void * data;
} utree_t;

typedef struct rtree_s
{
  char * label;
  double length;
  struct rtree_s * left;
  struct rtree_s * right;
  struct rtree_s * parent;
  int leaves;

  void * data;
} rtree_t;

typedef struct spec_array_entry
{
  double sum_speciation_edges_subtree;
  double coalescent_value;
  double speciation_value;
  double score_multi;
  double score_single;
  int taken_left_index;
  int taken_right_index;
} spec_entry;

typedef struct node_information_ptpmulti
{
  int num_edges_subtree;
  double sum_edges_subtree;
  double coalescent;
  spec_entry * spec_array;

  // additional data
  int num_known_speciation_edges;
  double sum_known_speciation_edges;
} node_information;

typedef struct node_information_score
{
  bool marked; // is the node a most recent common ancestor (mrca) node or not?
  int current_species_real; // for finding the "real" mrca
  int current_species_input; // for finding the alternative mrca
} score_information;

/* macros */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* options */

extern int opt_quiet;
extern char * opt_treefile;
extern char * opt_outfile;
extern char * opt_outgroup;
extern long opt_help;
extern long opt_version;
extern long opt_treeshow;

/* common data */

extern char errmsg[200];

extern long mmx_present;
extern long sse_present;
extern long sse2_present;
extern long sse3_present;
extern long ssse3_present;
extern long sse41_present;
extern long sse42_present;
extern long popcnt_present;
extern long avx_present;
extern long avx2_present;

/* functions in util.c */

void fatal(const char * format, ...);
void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned int progress);
void progress_done();
void * xmalloc(size_t size);
void * xrealloc(void *ptr, size_t size);
char * xstrchrnul(char *s, int c);
char * xstrdup(const char * s);
char * xstrndup(const char * s, size_t len);
long getusec(void);
void show_rusage();

/* functions in delimit.c */

void args_init(int argc, char ** argv);
void cmd_help();
void getentirecommandline(int argc, char * argv[]);
void fillheader();
void show_header();
void cmd_divtimes();

/* functions in parse_rtree.y */

rtree_t * rtree_parse_newick(const char * filename);

void rtree_destroy(rtree_t * root);

/* functions in parse_utree.y */

utree_t * utree_parse_newick(const char * filename,
                             int * tip_count);

void utree_destroy(utree_t * root);

/* functions in utree.c */

void utree_show_ascii(utree_t * tree);

char * utree_export_newick(utree_t * root);

int utree_traverse(utree_t * root,
                   int (*cbtrav)(utree_t *),
                   utree_t ** outbuffer);

int utree_query_tipnodes(utree_t * root,
                         utree_t ** node_list);

int utree_query_innernodes(utree_t * root,
                           utree_t ** node_list);

rtree_t * utree_convert_rtree(utree_t * root, int tip_count);

/* functions in rtree.c */

void rtree_show_ascii(rtree_t * tree);

char * rtree_export_newick(rtree_t * root);

int rtree_traverse(rtree_t * root,
                   int (*cbtrav)(rtree_t *),
                   rtree_t ** outbuffer);

int rtree_query_tipnodes(rtree_t * root,
                         rtree_t ** node_list);

int rtree_query_innernodes(rtree_t * root,
                           rtree_t ** node_list);

void rtree_reset_leaves(rtree_t * root);

/* functions in parse_rtree.y */

rtree_t * rtree_parse_newick(const char * filename);

/* functions in lca_utree.c */

void lca_init(utree_t * root);
utree_t * lca_compute(utree_t * tip1, utree_t * tip2);
void lca_destroy();

/* functions in arch.c */

unsigned long arch_get_memused();
unsigned long arch_get_memtotal();

/* functions in ptp_multi.c */

void ptp_multi_heuristic(rtree_t * rtree, bool multiple_lambda);

/* functions in score.c */

void compare_delimitation_tree(char * scorefile, rtree_t * tree);
