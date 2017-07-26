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

#define _GNU_SOURCE

#include <assert.h>
#include <ctype.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdbool.h>
#include <time.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if (defined(HAVE_CONFIG_H) && defined(HAVE_LIBGSL))
#include <gsl/gsl_cdf.h>
#endif

/* constants */

#define PROG_NAME PACKAGE
#define PROG_VERSION PACKAGE_VERSION

#ifdef __PPC__

#ifdef __LITTLE_ENDIAN__
#define PROG_CPU "ppc64le"
#else
#error "Big endian ppc64 CPUs not supported"
#endif

#else

#define PROG_CPU "x86_64"

#endif

#ifdef __APPLE__
#define PROG_OS "osx"
#include <sys/resource.h>
#include <sys/sysctl.h>
#endif

#ifdef __linux__
#define PROG_OS "linux"
#include <sys/resource.h>
#include <sys/sysinfo.h>
#endif

#ifdef _WIN32
#define PROG_OS "win"
#include <windows.h>
#include <psapi.h>
#endif

#define PROG_ARCH PROG_OS "_" PROG_CPU

#define PLL_FAILURE  0
#define PLL_SUCCESS  1
#define PLL_LINEALLOC 2048
#define PLL_ERROR_FILE_OPEN              1
#define PLL_ERROR_FILE_SEEK              2
#define PLL_ERROR_FILE_EOF               3
#define PLL_ERROR_FASTA_ILLEGALCHAR      4
#define PLL_ERROR_FASTA_UNPRINTABLECHAR  5
#define PLL_ERROR_FASTA_INVALIDHEADER    6
#define PLL_ERROR_MEM_ALLOC              7

#define LINEALLOC 2048

#define EVENT_SPECIATION 0
#define EVENT_COALESCENT 1

#define PTP_METHOD_SINGLE       0
#define PTP_METHOD_MULTI        1

#define REGEX_REAL   "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)"

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;

typedef struct dp_vector_s
{
  /* sum of speciation edge lengths of current subtree */
  double spec_edgelen_sum;

  /* coalescent logl of subtree for multi lambda */
  double coal_multi_logl;

  /* best single- and multi-rate log-likelihood for current subtree */
  double score_multi;
  double score_single;

  /* back-tracking information */
  int vec_left;
  int vec_right;

  unsigned int species_count;
  int filled;
} dp_vector_t;

typedef struct utree_s
{
  char * label;
  double length;
  int height;
  struct utree_s * next;
  struct utree_s * back;

  void * data;

  /* for finding the lca */
  int mark;

} utree_t;

typedef struct rtree_s
{
  char * label;
  double length;
  struct rtree_s * left;
  struct rtree_s * right;
  struct rtree_s * parent;
  int leaves;

  /* number of edges within current subtree with lengths greater than opt_minbr
     and corresponding sum */
  int edge_count;
  double edgelen_sum;
  double coal_logl;

  /* minimum number of speciation edges if current node is the start of a
     coalescent event, and the respective sum of lengths  */
  int spec_edge_count;
  double spec_edgelen_sum;

  /* which process does this node belong to (coalesent or speciation) */
  int event;

  /* slot in which the node resides when doing mcmc analysis */
  long mcmc_slot;
  long speciation_start;
  long speciation_count;
  double aic_weight_start;
  double aic_support;
  double support;

  /* dynamic programming vector */
  dp_vector_t * vector;

  /* auxialiary data */
  void * data;

  /* for generating random delimitations */
  int max_species_count;

  /* mark */
  int mark;
  char * sequence;

} rtree_t;

typedef struct pll_fasta
{
  FILE * fp;
  char line[LINEALLOC];
  const unsigned int * chrstatus;
  long no;
  long filesize;
  long lineno;
  long stripped_count;
  long stripped[256];
} pll_fasta_t;

typedef struct list_item_s
{
  void * data;
  struct list_item_s * next;
} list_item_t;

typedef struct list_s
{
  list_item_t * head;
  list_item_t * tail;
  long count;
} list_t;

typedef struct ht_item_s
{
  unsigned long key;
  void * value;
} ht_item_t;

typedef struct hashtable_s
{
  unsigned long table_size;
  unsigned long entries_count;
  list_t ** entries;
} hashtable_t;

typedef struct pair_s
{
  char * label;
  size_t index;
} pair_t;

/* macros */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* options */

extern int opt_quiet;
extern int opt_precision;
extern int opt_svg_showlegend;
extern long opt_help;
extern long opt_version;
extern long opt_treeshow;
extern long opt_mcmc_sample;
extern long opt_mcmc_steps;
extern long opt_mcmc_log;
extern long opt_mcmc_startml;
extern long opt_mcmc_startnull;
extern long opt_mcmc_startrandom;
extern long opt_mcmc_burnin;
extern long opt_mcmc_runs;
extern long opt_seed;
extern long opt_mcmc;
extern long opt_ml;
extern long opt_multi;
extern long opt_single;
extern long opt_method;
extern long opt_crop;
extern long opt_svg;
extern long opt_svg_width;
extern long opt_svg_fontsize;
extern long opt_svg_tipspace;
extern long opt_svg_marginleft;
extern long opt_svg_marginright;
extern long opt_svg_margintop;
extern long opt_svg_marginbottom;
extern long opt_svg_inner_radius;
extern double opt_mcmc_credible;
extern double opt_svg_legend_ratio;
extern double opt_pvalue;
extern double opt_minbr;
extern char * opt_treefile;
extern char * opt_outfile;
extern char * opt_outgroup;
extern char * opt_pdist_file;
extern char * cmdline;

/* common data */

extern char errmsg[200];

extern int pll_errno;
extern unsigned short global_xsubi[3];
extern const unsigned int pll_map_nt[256];
extern const unsigned int pll_map_fasta[256];

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

void fatal(const char * format, ...) __attribute__ ((noreturn));
void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned int progress);
void progress_done(void);
void * xmalloc(size_t size);
void * xcalloc(size_t nmemb, size_t size);
void * xrealloc(void *ptr, size_t size);
char * xstrchrnul(char *s, int c);
char * xstrdup(const char * s);
char * xstrndup(const char * s, size_t len);
long getusec(void);
FILE * xopen(const char * filename, const char * mode);
void random_init(unsigned short * rstate, long seedval);
double mptp_erand48(unsigned short * rstate);
long mptp_nrand48(unsigned short * rstate); 

/* functions in mptp.c */

void args_init(int argc, char ** argv);
void cmd_help(void);
void getentirecommandline(int argc, char * argv[]);
void fillheader(void);
void show_header(void);
void cmd_ml(void);
void cmd_multirun(void);
void cmd_auto(void);

/* functions in parse_rtree.y */

rtree_t * rtree_parse_newick(const char * filename);
void rtree_destroy(rtree_t * root);

/* functions in parse_utree.y */

utree_t * utree_parse_newick(const char * filename, unsigned int * tip_count);

void utree_destroy(utree_t * root);

/* functions in utree.c */

void utree_show_ascii(utree_t * tree);
char * utree_export_newick(utree_t * root);
int utree_query_tipnodes(utree_t * root, utree_t ** node_list);
int utree_query_innernodes(utree_t * root, utree_t ** node_list);
rtree_t * utree_convert_rtree(utree_t * root);
int utree_traverse(utree_t * root,
                   int (*cbtrav)(utree_t *),
                   utree_t ** outbuffer);
utree_t * utree_longest_branchtip(utree_t * node, unsigned int tip_count);
utree_t * utree_outgroup_lca(utree_t * root, unsigned int tip_count);
rtree_t * utree_crop(utree_t * lca);

/* functions in rtree.c */

void rtree_show_ascii(rtree_t * tree);
char * rtree_export_newick(rtree_t * root);
int rtree_query_tipnodes(rtree_t * root, rtree_t ** node_list);
int rtree_query_innernodes(rtree_t * root, rtree_t ** node_list);
void rtree_reset_info(rtree_t * root);
void rtree_print_tips(rtree_t * node, FILE * out);
int rtree_traverse(rtree_t * root,
                   int (*cbtrav)(rtree_t *),
                   unsigned short * rstate,
                   rtree_t ** outbuffer);
rtree_t * rtree_clone(rtree_t * node, rtree_t * parent);
int rtree_traverse_postorder(rtree_t * root,
                             int (*cbtrav)(rtree_t *),
                             rtree_t ** outbuffer);
rtree_t * get_outgroup_lca(rtree_t * root);
rtree_t * rtree_lca(rtree_t * root,
                    rtree_t ** tip_nodes,
                    unsigned int count);
rtree_t * rtree_crop(rtree_t * root, rtree_t * crop_root);
int rtree_height(rtree_t * root);

/* functions in parse_rtree.y */

rtree_t * rtree_parse_newick(const char * filename);

/* functions in lca_utree.c */

void lca_init(utree_t * root);
utree_t * lca_compute(utree_t * tip1, utree_t * tip2);
void lca_destroy(void);

/* functions in arch.c */

unsigned long arch_get_memused(void);
unsigned long arch_get_memtotal(void);
long arch_get_cores(void);

/* functions in dp.c */

void dp_init(rtree_t * tree);
void dp_free(rtree_t * tree);
void dp_ptp(rtree_t * rtree, long method);
void dp_set_pernode_spec_edges(rtree_t * node);

/* functions in svg.c */

void cmd_svg(rtree_t * rtree, long seed, const char * ext);

/* functions in likelihood.c */

double loglikelihood(long edge_count, double edgelen_sum);
int lrt(double nullmodel_logl, double ptp_logl, unsigned int df, double * pvalue);
double aic(double logl, long k, long n);

/* functions in output.c */

void output_info(FILE * out,
  		 long method,
		 double nullmodel_logl,
		 double logl,
		 double pvalue,
		 int lrt_result,
                 rtree_t * root,
                 unsigned int species_count);

FILE * open_file_ext(const char * extension, long seed);

/* functions in svg_landscape.c */

void svg_landscape(double mcmc_min_log, double mcmc_max_logl, long seed);
void svg_landscape_combined(double mcmc_min_log, double mcmc_max_logl, long runs, long * seed);

/* functions in random.c */

double random_delimitation(rtree_t * root,
                           long * delimited_species,
                           long * coal_edge_count,
                           double * coal_edgelen_sum,
                           long * spec_edge_count,
                           double * spec_edgelen_sum,
                           double * coal_score,
                           unsigned short * rstate);

/* functions in multirun.c */

void multirun(rtree_t * root, long method);

/* functions in fasta.c */

pll_fasta_t * pll_fasta_open(const char * filename,
                                        const unsigned int * map);

int pll_fasta_getnext(pll_fasta_t * fd, char ** head,
                                 long * head_len,  char ** seq,
                                 long * seq_len, long * seqno);

void pll_fasta_close(pll_fasta_t * fd);

long pll_fasta_getfilesize(pll_fasta_t * fd);

long pll_fasta_getfilepos(pll_fasta_t * fd);

int pll_fasta_rewind(pll_fasta_t * fd);

/* functions in auto.c */

void detect_min_bl(rtree_t * rtree);

/* functions in aic.c */

void aic_mcmc(rtree_t * tree,
              long method,
              unsigned short * rstate,
              long seed,
              double * mcmc_min_logl,
              double * mcmc_max_logl);

/* functions in hash.c */

unsigned long hash_djb2a(char * s);

unsigned long hash_fnv(char * s);

int hashtable_strcmp(void * x, void * y);

int hashtable_ptrcmp(void * x, void * y);

int hashtable_paircmp(void * stored, void * query);

void * hashtable_find(hashtable_t * ht,
                      void * x,
                      unsigned long hash,
                      int (*cb_cmp)(void *, void *));

hashtable_t * hashtable_create(unsigned long items_count);

int hashtable_insert(hashtable_t * ht,
                     void * x,
                     unsigned long hash,
                     int (*cb_cmp)(void *, void *));


/* functions in list.c */

void list_append(list_t * list, void * data);

void list_prepend(list_t * list, void * data);

void list_clear(list_t * list, void (*cb_dealloc)(void *));

list_t * list_create(void * data);

void hashtable_destroy(hashtable_t * ht, void (*cb_dealloc)(void *));
