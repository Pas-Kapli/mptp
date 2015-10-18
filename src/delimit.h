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
#include <regex.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdbool.h>
#include <gsl/gsl_cdf.h>

/* constants */

#define PROG_NAME "delimit"
#define PROG_VERSION "v0.0.1"

#ifdef __APPLE__
#define PROG_ARCH "macosx_x86_64"
#else
#define PROG_ARCH "linux_x86_64"
#endif

#define EVENT_SPECIATION 0
#define EVENT_COALESCENT 1

#define PRIOR_NONE              0
#define PRIOR_UNI               1
#define PRIOR_NBIN              2
#define PRIOR_BINOMIAL          3
#define PRIOR_GAMMA             4
#define PRIOR_DIRICHLET         5
#define PRIOR_BETA              6
#define PRIOR_EXP               7

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

  /* dynamic programming vector */
  dp_vector_t * vector;

  /* auxialiary data */
  void * data;
} rtree_t;


/* priors and hyperpriors */

typedef struct prior_s
{
  int dist;
  void * params;
} prior_t;

typedef struct gamma_params_s
{
  double alpha;
  double beta;
} gamma_params_t;

typedef struct beta_params_s
{
  double alpha;
  double beta;
} beta_params_t;

typedef struct bin_params_s
{
  int trials;
  double prob;
} bin_params_t;

typedef struct nbin_params_s
{
  double prob;
  int failures;
} nbin_params_t;

typedef struct dir_params_s 
{
  int k;
  unsigned int * a;
} dir_params_t;

typedef struct uni_params_s
{
  double min;
  double max;
} uni_params_t;

typedef struct exp_params_s
{
  double rate;
} exp_params_t;

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
extern long opt_ml_multi;
extern long opt_ml_single;
extern long opt_bayes_multi;
extern long opt_bayes_single;
extern long opt_ks_single;
extern long opt_ks_multi;
extern long opt_bayes_runs;
extern long opt_svg;
extern long opt_svg_width;
extern long opt_svg_fontsize;
extern long opt_svg_tipspace;
extern long opt_svg_marginleft;
extern long opt_svg_marginright;
extern long opt_svg_margintop;
extern long opt_svg_marginbottom;
extern long opt_svg_inner_radius;
extern double opt_svg_legend_ratio;
extern double opt_pvalue;
extern double opt_minbr;
extern char * opt_treefile;
extern char * opt_outfile;
extern char * opt_outgroup;
extern char * opt_scorefile;
extern char * cmdline;
extern prior_t * opt_prior;

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
int extract2f(char * text, double * a, double * b);
FILE * xopen(const char * filename, const char * mode);

/* functions in delimit.c */

void args_init(int argc, char ** argv);
void cmd_help(void);
void getentirecommandline(int argc, char * argv[]);
void fillheader(void);
void show_header(void);
void cmd_ml_single(void);
void cmd_ml_multi(void);
void cmd_score(void);
void cmd_ks_multi(void);
void cmd_ks_single(void);

/* functions in parse_rtree.y */

rtree_t * rtree_parse_newick(const char * filename);
void rtree_destroy(rtree_t * root);

/* functions in parse_utree.y */

utree_t * utree_parse_newick(const char * filename, int * tip_count);

void utree_destroy(utree_t * root);

/* functions in utree.c */

void utree_show_ascii(utree_t * tree);
char * utree_export_newick(utree_t * root);
int utree_query_tipnodes(utree_t * root, utree_t ** node_list);
int utree_query_innernodes(utree_t * root, utree_t ** node_list);
rtree_t * utree_convert_rtree(utree_t * root, int tip_count);
int utree_traverse(utree_t * root,
                   int (*cbtrav)(utree_t *),
                   utree_t ** outbuffer);

/* functions in rtree.c */

void rtree_show_ascii(rtree_t * tree);
char * rtree_export_newick(rtree_t * root);
int rtree_query_tipnodes(rtree_t * root, rtree_t ** node_list);
int rtree_query_innernodes(rtree_t * root, rtree_t ** node_list);
void rtree_reset_info(rtree_t * root);
void rtree_print_tips(rtree_t * node, FILE * out);
int rtree_traverse(rtree_t * root,
                   int (*cbtrav)(rtree_t *),
                   rtree_t ** outbuffer);

/* functions in parse_rtree.y */

rtree_t * rtree_parse_newick(const char * filename);

/* functions in lca_utree.c */

void lca_init(utree_t * root);
utree_t * lca_compute(utree_t * tip1, utree_t * tip2);
void lca_destroy(void);

/* functions in arch.c */

unsigned long arch_get_memused(void);
unsigned long arch_get_memtotal(void);

/* functions in dp.c */

void dp_init(rtree_t * tree);
void dp_free(rtree_t * tree);
void dp_ptp(rtree_t * rtree, int multi, prior_t * prior);
void dp_set_pernode_spec_edges(rtree_t * node);

/* functions in score.c */

void score_delimitation_tree(char * scorefile, rtree_t * tree);

/* functions in svg.c */

void cmd_svg(rtree_t * rtree);

/* functions in priors.c */

double prior_score(unsigned int species_count, prior_t * prior);
void init_logn_table(unsigned int n);
double dir_logpdf(double * x, dir_params_t * params);
double gamma_logpdf(double x, gamma_params_t * params);
double beta_logpdf(double x, beta_params_t * params);
double bin_logpmf(unsigned int k, bin_params_t * params);
double nbin_logpmf(unsigned int k, nbin_params_t * params);
double uni_logpdf(double x, uni_params_t * params);
double exp_logpdf(double x, exp_params_t * params);

/* functions in knapsack.c */

void dp_knapsack(rtree_t * root, int method);

/* functions in likelihood.c */

double loglikelihood(int edge_count, double edgelen_sum);
int lrt(double nullmodel_logl, double ptp_logl, unsigned int df, double * pvalue);

/* functions in output.c */

void output_info(FILE * out,
  		 int method,
		 double nullmodel_logl,
		 double logl,
		 double pvalue,
		 int lrt_result,
                 rtree_t * root,
                 int species_count);

FILE * open_file_ext(const char * extension);
