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

#include "mptp.h"

static char * progname;
static char progheader[80];
char * cmdline;

/* global error message buffer */
char errmsg[200] = {0};

/* number of mandatory options for the user to input */
static const char mandatory_options_count = 2;
static const char * mandatory_options_list = " --tree_file --output_file";

/* options */
int pll_errno;
int opt_quiet;
int opt_precision;
int opt_svg_showlegend;
long opt_help;
long opt_version;
long opt_treeshow;
long opt_ml_multi;
long opt_ml_single;
long opt_bayes_multi;
long opt_bayes_single;
long opt_bayes_sample;
long opt_bayes_runs;
long opt_bayes_log;
long opt_bayes_startnull;
long opt_bayes_startrandom;
long opt_bayes_burnin;
long opt_bayes_chains;
long opt_seed;
long opt_crop;
long opt_svg;
long opt_svg_width;
long opt_svg_fontsize;
long opt_svg_tipspace;
long opt_svg_marginleft;
long opt_svg_marginright;
long opt_svg_margintop;
long opt_svg_marginbottom;
long opt_svg_inner_radius;
double opt_svg_legend_ratio;
double opt_pvalue;
double opt_minbr;
char * opt_treefile;
char * opt_outfile;
char * opt_outgroup;
char * opt_scorefile;
char * opt_pdist_file;
prior_t * opt_prior;

static void dealloc_prior()
{
  if (opt_prior) 
  {
    if (opt_prior->params)
      free(opt_prior->params);
    free(opt_prior);
  }
}

static struct option long_options[] =
{
  {"help",               no_argument,       0, 0 },  /*  0 */
  {"version",            no_argument,       0, 0 },  /*  1 */
  {"quiet",              no_argument,       0, 0 },  /*  2 */
  {"tree_file",          required_argument, 0, 0 },  /*  3 */
  {"tree_show",          no_argument,       0, 0 },  /*  4 */
  {"output_file",        required_argument, 0, 0 },  /*  5 */
  {"ml_multi",           no_argument,       0, 0 },  /*  6 */
  {"ml_single",          no_argument,       0, 0 },  /*  7 */
  {"outgroup",           required_argument, 0, 0 },  /*  8 */
  {"score",              required_argument, 0, 0 },  /*  9 */
  {"pvalue",             required_argument, 0, 0 },  /* 10 */
  {"minbr",              required_argument, 0, 0 },  /* 11 */
  {"svg_width",          required_argument, 0, 0 },  /* 12 */
  {"svg_fontsize",       required_argument, 0, 0 },  /* 13 */
  {"svg_tipspacing",     required_argument, 0, 0 },  /* 14 */
  {"svg_legend_ratio",   required_argument, 0, 0 },  /* 15 */
  {"svg_nolegend",       no_argument,       0, 0 },  /* 16 */
  {"svg_marginleft",     required_argument, 0, 0 },  /* 17 */
  {"svg_marginright",    required_argument, 0, 0 },  /* 18 */
  {"svg_margintop",      required_argument, 0, 0 },  /* 19 */
  {"svg_marginbottom",   required_argument, 0, 0 },  /* 20 */
  {"svg_inner_radius",   required_argument, 0, 0 },  /* 21 */
  {"precision",          required_argument, 0, 0 },  /* 22 */
  {"bayes_multi",        required_argument, 0, 0 },  /* 23 */
  {"bayes_single",       required_argument, 0, 0 },  /* 24 */
  {"prior_exp",          required_argument, 0, 0 },  /* 25 */
  {"prior_gamma",        required_argument, 0, 0 },  /* 26 */
  {"bayes_sample",       required_argument, 0, 0 },  /* 27 */
  {"bayes_log",          no_argument,       0, 0 },  /* 28 */
  {"seed",               required_argument, 0, 0 },  /* 29 */
  {"bayes_startnull",    no_argument,       0, 0 },  /* 30 */ 
  {"bayes_burnin",       required_argument, 0, 0 },  /* 31 */
  {"bayes_startrandom",  no_argument,       0, 0 },  /* 32 */
  {"bayes_chains",       required_argument, 0, 0 },  /* 33 */
  {"minbr_auto",         required_argument, 0, 0 },  /* 34 */
  {"outgroup_crop",      no_argument,       0, 0 },  /* 35 */
  { 0, 0, 0, 0 }
};

void args_init(int argc, char ** argv)
{
  int option_index = 0;
  int c;
  int mand_options = 0;

  /* set defaults */

  progname = argv[0];

  opt_help = 0;
  opt_version = 0;
  opt_treeshow = 0;
  opt_treefile = NULL;
  opt_outfile = NULL;
  opt_outgroup = NULL;
  opt_pdist_file = NULL;
  opt_quiet = 0;
  opt_ml_multi = 0;
  opt_ml_single = 0;
  opt_bayes_multi = 0;
  opt_bayes_single = 0;
  opt_scorefile = NULL;
  opt_pvalue = 0.001;
  opt_minbr = 0.0001;
  opt_precision = 7;
  opt_bayes_runs = 1;
  opt_bayes_sample = 1000;
  opt_bayes_startnull = 0;
  opt_bayes_startrandom = 0;
  opt_bayes_log = 0;
  opt_bayes_burnin = 1;
  opt_bayes_chains = 0;
  opt_seed = (long)time(NULL);
  opt_crop = 0;

  opt_svg_width = 1920;
  opt_svg_fontsize = 12;
  opt_svg_tipspace = 20;
  opt_svg_legend_ratio = 0.1;
  opt_svg_showlegend = 1;
  opt_svg_marginleft = 20;
  opt_svg_marginright = 20;
  opt_svg_margintop = 20;
  opt_svg_marginbottom = 20;
  opt_svg_inner_radius = 0;

  opt_prior = NULL;

  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) == 0)
  {
    char * end;
    switch (option_index)
    {
      case 0:
        opt_help = 1;
        break;

      case 1:
        opt_version = 1;
        break;

      case 2:
        opt_quiet = 1;
        break;

      case 3:
        free(opt_treefile);
        opt_treefile = optarg;
        break;

      case 4:
        opt_treeshow = 1;
        break;

      case 5:
        opt_outfile = optarg;
        break;

      case 6:
        opt_ml_multi = 1;
        break;

      case 7:
        opt_ml_single = 1;
        break;

      case 8:
        opt_outgroup = optarg;
        break;

      case 9:
        free(opt_scorefile);
        opt_scorefile = optarg;
        break;

      case 10:
        opt_pvalue = strtod(optarg, &end);
        if (end == optarg) {
          fatal(" is not a valid number.\n");
        }
        break;

      case 11:
        opt_minbr = strtod(optarg, &end);
        if (end == optarg) {
          fatal(" is not a valid number.\n");
        }
        break;

      case 12:
        opt_svg_width = atoi(optarg);
        break;

      case 13:
        opt_svg_fontsize = atol(optarg);
        break;

      case 14:
        opt_svg_tipspace = atol(optarg);
        break;

      case 15:
        opt_svg_legend_ratio = atof(optarg);
        break;

      case 16:
        opt_svg_showlegend = 0;
        break;

      case 17:
        opt_svg_marginleft = atol(optarg);
        break;

      case 18:
        opt_svg_marginright = atol(optarg);
        break;

      case 19:
        opt_svg_margintop = atol(optarg);
        break;

      case 20:
        opt_svg_marginbottom = atol(optarg);
        break;

      case 21:
        opt_svg_inner_radius = atol(optarg);
        break;

      case 22:
        opt_precision = atoi(optarg);
        break;

      case 23:
        opt_bayes_multi = 1;
        opt_bayes_runs = atoi(optarg);
        break;

      case 24:
        opt_bayes_single = 1;
        opt_bayes_runs = atol(optarg);
        break;

      case 25:
        dealloc_prior();

        opt_prior = (prior_t *)xmalloc(sizeof(prior_t));
        opt_prior->dist = PRIOR_EXP;

        exp_params_t * exp_params =(exp_params_t *)xmalloc(sizeof(exp_params_t));
        exp_params->rate = atof(optarg);
        opt_prior->params = (void *)exp_params;
        break;

      case 26:
        dealloc_prior();

        opt_prior = (prior_t *)xmalloc(sizeof(prior_t));
        opt_prior->dist = PRIOR_GAMMA;

        gamma_params_t * gamma_params = (gamma_params_t *)xmalloc(sizeof(gamma_params_t));
        if (!extract2f(optarg, &(gamma_params->alpha), &(gamma_params->beta)))
          fatal("Incorrect format for --prior_gamma");
        opt_prior->params = (void *)gamma_params;
        break;

      case 27:
        opt_bayes_sample = atol(optarg);
        break;

      case 28:
        opt_bayes_log = 1;
        break;

      case 29:
        opt_seed = atol(optarg);
        break;

      case 30:
        opt_bayes_startnull = 1;
        break;

      case 31:
        opt_bayes_burnin = atol(optarg);
        break;

      case 32:
        opt_bayes_startrandom = 1;
        break;

      case 33:
        opt_bayes_chains = atol(optarg);
        break;

      case 34:
        free(opt_pdist_file);
        opt_pdist_file = optarg;
        break;

      case 35:
        opt_crop = 1;
        break;

      default:
        fatal("Internal error in option parsing");
    }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands  = 0;

  /* check for mandatory options */
  if (opt_treefile)
    mand_options++;
  if (opt_outfile)
    mand_options++;

  /* check for number of independent commands selected */
  if (opt_version)
    commands++;
  if (opt_help)
    commands++;
  if (opt_ml_multi)
    commands++;
  if (opt_ml_single)
    commands++;
  if (opt_bayes_multi)
    commands++;
  if (opt_bayes_single)
    commands++;
  if (opt_scorefile)
    commands++;
  if (opt_pdist_file)
    commands++;

  /* if more than one independent command, fail */
  if (commands > 1)
    fatal("More than one command specified");

  /* if no command specified, turn on --help */
  if (!commands)
  {
    opt_help = 1;
    return;
  }
  /* check for mandatory options */
  if (!opt_version && !opt_help)
    if (mand_options != mandatory_options_count)
      fatal("Mandatory options are:\n\n%s", mandatory_options_list);

}

void cmd_help()
{
  fprintf(stderr,
          "Usage: %s [OPTIONS]\n", progname);
  fprintf(stderr,
          "\n"
          "General options:\n"
          "  --help                    display help information.\n"
          "  --version                 display version information.\n"
          "  --tree_show               display an ASCII version of the tree.\n"
          "  --ml_multi                Maximum-likelihood with one lambda per coalescent.\n"
          "  --ml_single               Maximum-likelihood with one lambda for all coalescent.\n"
          "  --bayes_multi INT         Bayesian with own lambda per coalescent (INT runs).\n"
          "  --bayes_single INT        Bayesian with one lambda for all coalescent (INT steps).\n"
          "  --bayes_sample INT        Sample every INT iteration (default: 1000).\n"
          "  --bayes_log               Log samples and create SVG plot of log-likelihoods.\n"
          "  --bayes_burnin INT        Ignore all MCMC steps below threshold.\n"
          "  --bayes_chains INT        Run multiple chains.\n"
          "  --score                   Compare given species delimitation with optimal one induced by the tree.\n"
          "  --pvalue                  Set p-value for LRT (default: 0.001)\n"
          "  --minbr REAL              Set minimum branch length (default: 0.0001)\n"
          "  --minbr_auto FILENAME     Detect minimum branch length from FASTA p-distances\n"
          "  --outgroup TAXON          If input tree is unrooted, use TAXON as outgroup (default: taxon with longest branch).\n"
          "  --quiet                   only output warnings and fatal errors to stderr.\n"
          "  --precision               Precision of decimal part of floating point numbers on output (default: 7).\n"
          "  --seed                    Seed for pseudo-random number generator.\n"
          "Prior options:\n"
          "  --prior_exp REAL          Rate of exponential prior.\n"
          "  --prior_ln REAL,REAL      Log-normal prior with mean (first param) and standard deviation (second param).\n"
          "  --prior_uni REAL,REAL     Uniform prior with minimum (first param) and maximum (second param) bounds.\n"
          "  --prior_bin INT,REAL      Binomial prior with number of trials (first param) and success probability (second param).\n"
          "  --prior_nbin INT,REAL     Negative binomial prior with number of failures (first param) and success probability (second param).\n"
          "  --prior_gamma REAL,REAL   Gamma distribution with shape (first param) and rate (second param).\n"
          "  --prior_beta REAL,REAL    Beta distribution with alpha shape (first param) and beta shape (second param).\n"
          "Input and output options:\n"
          "  --tree_file FILENAME      tree file in newick format.\n"
          "  --output_file FILENAME    output file name.\n"
          "Visualization options:\n"
          "  --svg_width INT           Width of SVG tree in pixels (default: 1920).\n"
          "  --svg_fontsize INT        Size of font in SVG image. (default: 12)\n"
          "  --svg_tipspacing INT      Vertical space between taxa in SVG tree (default: 20).\n"
          "  --svg_legend_ratio <0..1> Ratio of total tree length to be displayed as legend line.\n"
          "  --svg_nolegend            Hides legend.\n"
          "  --svg_marginleft          Left margin in pixels (default: 20).\n"
          "  --svg_marginright         Right margin in pixels (default: 20).\n"
          "  --svg_margintop           Top margin in pixels (default: 20).\n"
          "  --svg_marginbottom        Bottom margin in pixels (default: 20).\n"
          "  --svg_inner_radius        Radius of inner nodes in pixels (default: 0).\n"
         );
}

static rtree_t * load_tree(void)
{
  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  rtree_t * rtree = rtree_parse_newick(opt_treefile);

  if (!rtree)
  {
    int tip_count;
    utree_t * utree = utree_parse_newick(opt_treefile, &tip_count);
    if (!utree)
      fatal("Tree is neither unrooted nor rooted. Go fix your tree.");

    if (!opt_quiet)
    {
      fprintf(stdout, "Loaded unrooted tree...\n");
      fprintf(stdout, "Converting to rooted tree...\n");
    }

    /* if outgroup was not specified, get the node with the longest branch */
    utree_t * og_root = NULL;

    /* if outgroup was not specified, get the tip with the longest branch */
    if (!opt_outgroup)
    {
      og_root = utree_longest_branchtip(utree, tip_count);
      assert(og_root);
      fprintf(stdout, 
              "Selected %s as outgroup based on longest tip-branch criterion\n",
              og_root->label);
    }
    else
    {
      /* get LCA of out group */
      og_root = utree_outgroup_lca(utree, tip_count);
      if (!og_root)
      {
        utree_destroy(utree);
        fatal("Outgroup must be a single tip or a list of all tips of a subtree");
      }
    }

    
    if (opt_crop)
    {
      rtree = utree_crop(og_root);
    }
    else
    {
      rtree = utree_convert_rtree(og_root);
    }

    utree_destroy(utree);
  }
  else
  {
    if (!opt_quiet)
      fprintf(stdout, "Loaded rooted tree...\n");

    if (opt_crop)
    {
      if (!opt_outgroup)
        fatal("--outgroup must be specified when using --outgroup_crop.");

      /* get LCA of outgroup */
      rtree_t * og_root = get_outgroup_lca(rtree);

      /* crop outgroup from tree */
      rtree = rtree_crop(rtree,og_root);
      if (!rtree)
        fatal("Cropping the outgroup leads to less than two tips.");
    }
  }

  return rtree;
}

void cmd_auto()
{
  rtree_t * rtree = load_tree();

  detect_min_bl(rtree);

  /* deallocate tree structure */
  rtree_destroy(rtree);
}

void cmd_ml_multi()
{

  rtree_t * rtree = load_tree();

  dp_init(rtree);
  dp_set_pernode_spec_edges(rtree);
  dp_ptp(rtree, PTP_METHOD_MULTI, opt_prior);
  dp_free(rtree);

  if (opt_treeshow)
    rtree_show_ascii(rtree);

  cmd_svg(rtree, opt_seed);

  /* deallocate tree structure */
  rtree_destroy(rtree);

  if (!opt_quiet)
    fprintf(stdout, "Done...\n");
}

void cmd_multichain(int method)
{
  rtree_t * rtree = load_tree();

  /* init random number generator */
  srand48(opt_seed);

  multichain(rtree, method, opt_prior);

  if (opt_treeshow)
    rtree_show_ascii(rtree);
   
  if (!opt_quiet)
    fprintf(stdout, "Done...\n");

}

void cmd_ml_single()
{

  rtree_t * rtree = load_tree();

  dp_init(rtree);
  dp_set_pernode_spec_edges(rtree);
  dp_ptp(rtree, PTP_METHOD_SINGLE, opt_prior);
  dp_free(rtree);

  if (opt_treeshow)
    rtree_show_ascii(rtree);

  cmd_svg(rtree, opt_seed);

  /* deallocate tree structure */
  rtree_destroy(rtree);

  if (!opt_quiet)
    fprintf(stdout, "Done...\n");
}

void cmd_bayes(int method)
{
  struct drand48_data rstate;
  double bayes_min_logl = 0;
  double bayes_max_logl = 0;

  if (opt_bayes_burnin < 1 || opt_bayes_burnin > opt_bayes_runs)
    fatal("--opt_bayes_burnin must be a positive integer smaller or equal to --opt_bayes_runs");

  if (opt_bayes_chains)
  {
    cmd_multichain(method);
    return;
  }


  rtree_t * rtree = load_tree();

  /* init random number generator */
  srand48_r(opt_seed, &rstate);

  dp_init(rtree);
  dp_set_pernode_spec_edges(rtree);
  bayes(rtree,
        method,
        opt_prior,
        &rstate,
        opt_seed,
        &bayes_min_logl,
        &bayes_max_logl);
  dp_free(rtree);

  if (opt_bayes_log)
    svg_landscape(bayes_min_logl, bayes_max_logl, opt_seed);

  if (opt_treeshow)
    rtree_show_ascii(rtree);
   
  char * newick = rtree_export_newick(rtree);

  if (!opt_quiet)
    fprintf(stdout,
            "Creating tree with support values in %s.%ld.tree ...\n",
            opt_outfile,
            opt_seed);

  FILE * newick_fp = open_file_ext("tree", opt_seed);
  fprintf(newick_fp, "%s\n", newick);
  fclose(newick_fp);

  cmd_svg(rtree, opt_seed);
  /* deallocate tree structure */
  rtree_destroy(rtree);

  free(newick);

  if (!opt_quiet)
    fprintf(stdout, "Done...\n");
}

void cmd_score()
{

  rtree_t * rtree = load_tree();

  /* TODO: Sarah's score function should be called here */
  score_delimitation_tree(opt_scorefile, rtree);

  if (opt_treeshow)
    rtree_show_ascii(rtree);

  cmd_svg(rtree, opt_seed);

  /* deallocate tree structure */
  rtree_destroy(rtree);

  if (!opt_quiet)
    fprintf(stdout, "Done...\n");
}

void getentirecommandline(int argc, char * argv[])
{
  int len = 0;
  int i;

  for (i = 0; i < argc; ++i)
    len += strlen(argv[i]);

  cmdline = (char *)xmalloc((size_t)(len + argc + 1));
  cmdline[0] = 0;

  for (i = 0; i < argc; ++i)
  {
    strcat(cmdline, argv[i]);
    strcat(cmdline, " ");
  }
}

void fillheader()
{
  snprintf(progheader, 80,
           "%s %s_%s, %1.fGB RAM, %ld cores",
           PROG_NAME, PROG_VERSION, PROG_ARCH,
           arch_get_memtotal() / 1024.0 / 1024.0 / 1024.0,
           sysconf(_SC_NPROCESSORS_ONLN));
}

void show_header()
{
  fprintf(stdout, "%s\n", progheader);
  fprintf(stdout, "https://github.com/Pas-Kapli/delimit\n");
  fprintf(stdout,"\n");
}

int main (int argc, char * argv[])
{
  fillheader();
  getentirecommandline(argc, argv);

  args_init(argc, argv);

  show_header();

  srand48(opt_seed);

  if (opt_help)
  {
    cmd_help();
  }
  else if (opt_ml_multi)
  {
    cmd_ml_multi();
  }
  else if (opt_ml_single)
  {
    cmd_ml_single();
  }
  else if (opt_bayes_multi)
  {
    cmd_bayes(PTP_METHOD_MULTI);
  }
  else if (opt_bayes_single)
  {
    cmd_bayes(PTP_METHOD_SINGLE);
  }
  else if (opt_scorefile)
  {
    cmd_score();
  }
  else if (opt_pdist_file)
  {
    cmd_auto();
  }

  /* free prior information */
  dealloc_prior();

  free(cmdline);
  return (0);
}
