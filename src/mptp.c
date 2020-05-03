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

/* global pseudo-random number generator 48-bit state */
unsigned short global_xsubi[3];

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
long opt_method;
long opt_mcmc_sample;
long opt_mcmc_steps;
long opt_mcmc_log;
long opt_mcmc_startnull;
long opt_mcmc_startrandom;
long opt_mcmc_startml;
long opt_mcmc_burnin;
long opt_mcmc_runs;
long opt_seed;
long opt_mcmc;
long opt_ml;
long opt_multi;
long opt_single;
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
double opt_mcmc_credible;
double opt_svg_legend_ratio;
double opt_pvalue;
double opt_minbr;
char * opt_treefile;
char * opt_outfile;
char * opt_outgroup;
char * opt_pdist_file;

static struct option long_options[] =
{
  {"help",               no_argument,       0, 0 },  /*  0 */
  {"version",            no_argument,       0, 0 },  /*  1 */
  {"quiet",              no_argument,       0, 0 },  /*  2 */
  {"tree_file",          required_argument, 0, 0 },  /*  3 */
  {"tree_show",          no_argument,       0, 0 },  /*  4 */
  {"output_file",        required_argument, 0, 0 },  /*  5 */
  {"outgroup",           required_argument, 0, 0 },  /*  6 */
  {"pvalue",             required_argument, 0, 0 },  /*  7 */
  {"minbr",              required_argument, 0, 0 },  /*  8 */
  {"svg_width",          required_argument, 0, 0 },  /*  9 */
  {"svg_fontsize",       required_argument, 0, 0 },  /* 10 */
  {"svg_tipspacing",     required_argument, 0, 0 },  /* 11 */
  {"svg_legend_ratio",   required_argument, 0, 0 },  /* 12 */
  {"svg_nolegend",       no_argument,       0, 0 },  /* 13 */
  {"svg_marginleft",     required_argument, 0, 0 },  /* 14 */
  {"svg_marginright",    required_argument, 0, 0 },  /* 15 */
  {"svg_margintop",      required_argument, 0, 0 },  /* 16 */
  {"svg_marginbottom",   required_argument, 0, 0 },  /* 17 */
  {"svg_inner_radius",   required_argument, 0, 0 },  /* 18 */
  {"precision",          required_argument, 0, 0 },  /* 19 */
  {"mcmc_sample",        required_argument, 0, 0 },  /* 20 */
  {"mcmc_log",           no_argument,       0, 0 },  /* 21 */
  {"seed",               required_argument, 0, 0 },  /* 22 */
  {"mcmc_startnull",     no_argument,       0, 0 },  /* 23 */
  {"mcmc_burnin",        required_argument, 0, 0 },  /* 24 */
  {"mcmc_startrandom",   no_argument,       0, 0 },  /* 25 */
  {"mcmc_runs",          required_argument, 0, 0 },  /* 26 */
  {"minbr_auto",         required_argument, 0, 0 },  /* 27 */
  {"outgroup_crop",      no_argument,       0, 0 },  /* 28 */
  {"mcmc_credible",      required_argument, 0, 0 },  /* 29 */
  {"mcmc",               required_argument, 0, 0 },  /* 30 */
  {"ml",                 no_argument,       0, 0 },  /* 31 */
  {"single",             no_argument,       0, 0 },  /* 32 */
  {"multi",              no_argument,       0, 0 },  /* 33 */
  {"mcmc_startml",       no_argument,       0, 0 },  /* 34 */
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
  opt_pvalue = 0.001;
  opt_minbr = 0.0001;
  opt_precision = 7;
  opt_mcmc_steps = 0;
  opt_mcmc_sample = 1000;
  opt_mcmc_startnull = 0;
  opt_mcmc_startrandom = 0;
  opt_mcmc_startml = 0;
  opt_mcmc_log = 0;
  opt_mcmc_burnin = 1;
  opt_mcmc_runs = 1;
  opt_mcmc_credible = 0.95;
  opt_seed = (long)time(NULL);
  opt_crop = 0;
  opt_ml = 0;
  opt_mcmc = 0;
  opt_method = PTP_METHOD_MULTI;
  opt_multi = 0;
  opt_single = 0;

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
        opt_outgroup = optarg;
        break;

      case 7:
        opt_pvalue = strtod(optarg, &end);
        if (end == optarg) {
          fatal(" is not a valid number.\n");
        }
        break;

      case 8:
        opt_minbr = strtod(optarg, &end);
        if (end == optarg) {
          fatal(" is not a valid number.\n");
        }
        break;

      case 9:
        opt_svg_width = atoi(optarg);
        break;

      case 10:
        opt_svg_fontsize = atol(optarg);
        break;

      case 11:
        opt_svg_tipspace = atol(optarg);
        break;

      case 12:
        opt_svg_legend_ratio = atof(optarg);
        break;

      case 13:
        opt_svg_showlegend = 0;
        break;

      case 14:
        opt_svg_marginleft = atol(optarg);
        break;

      case 15:
        opt_svg_marginright = atol(optarg);
        break;

      case 16:
        opt_svg_margintop = atol(optarg);
        break;

      case 17:
        opt_svg_marginbottom = atol(optarg);
        break;

      case 18:
        opt_svg_inner_radius = atol(optarg);
        break;

      case 19:
        opt_precision = atoi(optarg);
        break;

      case 20:
        opt_mcmc_sample = atol(optarg);
        break;

      case 21:
        opt_mcmc_log = 1;
        break;

      case 22:
        opt_seed = atol(optarg);
        break;

      case 23:
        opt_mcmc_startnull = 1;
        break;

      case 24:
        opt_mcmc_burnin = atol(optarg);
        break;

      case 25:
        opt_mcmc_startrandom = 1;
        break;

      case 26:
        opt_mcmc_runs = atol(optarg);
        break;

      case 27:
        free(opt_pdist_file);
        opt_pdist_file = optarg;
        break;

      case 28:
        opt_crop = 1;
        break;

      case 29:
        opt_mcmc_credible = atof(optarg);
        break;

      case 30:
        opt_mcmc = 1;
        opt_mcmc_steps = atol(optarg);
        break;

      case 31:
        opt_ml = 1;
        break;

      case 32:
        opt_method = PTP_METHOD_SINGLE;
        opt_single = 1;
        break;

      case 33:
        opt_method = PTP_METHOD_MULTI;
        opt_multi = 1;
        break;

      case 34:
        opt_mcmc_startml = 1;
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
  if (opt_pdist_file)
    commands++;
  if (opt_mcmc)
    commands++;
  if (opt_ml)
    commands++;

  /* if more than one independent command, fail */
  if (commands > 1)
  {
    if (!(commands == 2 && opt_pdist_file && (opt_ml || opt_mcmc)))
      fatal("More than one command specified");
  }

  /* if more than one independent command, fail */
  if (opt_mcmc_startrandom + opt_mcmc_startnull + opt_mcmc_startml > 1)
    fatal("You can only select one out of --mcmc_startrandom, --mcmc_startnull, --mcmc_startml");

  /* if more than one independent command, fail */
  if (opt_multi && opt_single)
    fatal("You can either specify --multi or --single, but not both at once.");

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
          "Examples:\n"
          "  mptp --ml --multi --tree_file tree.newick --output_file output\n"
          "  mptp --mcmc 50000000 --multi --mcmc_sample 1000000 --mcmc_burnin 1000000 --tree_file tree.newick --output_file output\n\n"
          "General options:\n"
          "  --help                    display help information.\n"
          "  --version                 display version information.\n"
          "  --tree_show               display an ASCII version of the tree.\n"
          "  --multi                   Use one lambda per coalescent (this is default).\n"
          "  --single                  Use one lambda for all coalescent.\n"
          "  --ml                      Maximum-likelihood heuristic.\n"
          "  --mcmc INT                Support values for the delimitation (INT steps).\n"
          "  --mcmc_sample INT         Sample every INT iteration (default: 1000).\n"
          "  --mcmc_log                Log samples and create SVG plot of log-likelihoods.\n"
          "  --mcmc_burnin INT         Ignore all MCMC steps below threshold.\n"
          "  --mcmc_runs INT           Perform multiple MCMC runs.\n"
          "  --mcmc_credible <0..1>    Credible interval (default: 0.95).\n"
          "  --mcmc_startnull          Start each run with the null model (one single species).\n"
          "  --mcmc_startrandom        Start each run with a random delimitation.\n"
          "  --mcmc_startml            Start each run with the delimitation obtained by the Maximum-likelihood heuristic.\n"
          "  --pvalue REAL             Set p-value for LRT (default: 0.001)\n"
          "  --minbr REAL              Set minimum branch length (default: 0.0001)\n"
          "  --minbr_auto FILENAME     Detect minimum branch length from FASTA p-distances\n"
          "  --outgroup TAXA           Root unrooted tree at outgroup (default: taxon with longest branch).\n"
          "  --outgroup_crop           Crop outgroup from tree\n"
          "  --quiet                   only output warnings and fatal errors to stderr.\n"
          "  --precision INT           Precision of floating point numbers on output (default: 7).\n"
          "  --seed                    Seed for pseudo-random number generator.\n"
          "\n"
          "Input and output options:\n"
          "  --tree_file FILENAME      tree file in newick format.\n"
          "  --output_file FILENAME    output file name.\n"
          "\n"
          "Visualization options:\n"
          "  --svg_width INT           Width of SVG tree in pixels (default: 1920).\n"
          "  --svg_fontsize INT        Size of font in SVG image. (default: 12)\n"
          "  --svg_tipspacing INT      Vertical space between taxa in SVG tree (default: 20).\n"
          "  --svg_legend_ratio <0..1> Ratio of total tree length to be displayed as legend line.\n"
          "  --svg_nolegend            Hides legend.\n"
          "  --svg_marginleft INT      Left margin in pixels (default: 20).\n"
          "  --svg_marginright INT     Right margin in pixels (default: 20).\n"
          "  --svg_margintop INT       Top margin in pixels (default: 20).\n"
          "  --svg_marginbottom  INT   Bottom margin in pixels (default: 20).\n"
          "  --svg_inner_radius INT    Radius of inner nodes in pixels (default: 0).\n"
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
    unsigned int tip_count;
    utree_t * utree = utree_parse_newick(opt_treefile, &tip_count);
    if (!utree)
      fatal("Tree is neither unrooted nor rooted.");

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

  detect_min_bl(rtree, 0);

  /* deallocate tree structure */
  rtree_destroy(rtree);
}

void cmd_ml(void)
{
  rtree_t * rtree = load_tree();
  if (opt_pdist_file)
    opt_minbr = detect_min_bl(rtree, 1);

  dp_init(rtree);
  dp_set_pernode_spec_edges(rtree);
  dp_ptp(rtree, opt_method);
  dp_free(rtree);

  if (opt_treeshow)
    rtree_show_ascii(rtree);

  cmd_svg(rtree, opt_seed, "svg");

  /* deallocate tree structure */
  rtree_destroy(rtree);

  if (!opt_quiet)
    fprintf(stdout, "Done...\n");
}

void cmd_multirun(void)
{
  if (opt_mcmc_steps == 0)
    fatal("The number of runs specified after --mcmc must be a positive integer greater than zero");

  if (opt_mcmc_burnin < 1 || opt_mcmc_burnin > opt_mcmc_steps)
    fatal("--opt_mcmc_burnin must be a positive integer smaller or equal to --opt_mcmc_steps");

  if (opt_mcmc_credible < 0 || opt_mcmc_credible > 1)
    fatal("--opt_mcmc_credible must be a real number between 0 and 1");

  rtree_t * rtree = load_tree();

  if (opt_pdist_file)
    opt_minbr = detect_min_bl(rtree, 1);

  multirun(rtree, opt_method);

  if (opt_treeshow)
    rtree_show_ascii(rtree);

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
           arch_get_cores());
}

void show_header()
{
  fprintf(stdout, "%s\n", progheader);
  fprintf(stdout, "https://github.com/Pas-Kapli/mptp\n");
  fprintf(stdout,"\n");
}

int main (int argc, char * argv[])
{
  fillheader();
  getentirecommandline(argc, argv);

  args_init(argc, argv);

  show_header();

  /* init random number generator and maintain compatibility with srand48 */
  random_init(global_xsubi,opt_seed);

  if (opt_help)
  {
    cmd_help();
  }
  else if (opt_pdist_file && !(opt_mcmc || opt_ml))
  {
    cmd_auto();
  }
  else if (opt_mcmc)
  {
    cmd_multirun();
  }
  else if (opt_ml)
  {
    cmd_ml();
  }

  free(cmdline);
  return (0);
}
