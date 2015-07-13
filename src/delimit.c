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

#include "delimit.h"

static char * progname;
static char progheader[80];
static char * cmdline;

/* global error message buffer */
char errmsg[200] = {0};

/* number of mandatory options for the user to input */
static const char mandatory_options_count = 2;
static const char * mandatory_options_list = " --tree_file --output_file";

/* options */
char * opt_treefile;
char * opt_outfile;
char * opt_outgroup;
int opt_quiet;
long opt_help;
long opt_version;
long opt_treeshow;
long opt_ptpmulti;
long opt_ptpsingle;


static struct option long_options[] =
{
  {"help",               no_argument,       0, 0 },  /*  0 */
  {"version",            no_argument,       0, 0 },  /*  1 */
  {"quiet",              no_argument,       0, 0 },  /*  2 */
  {"tree_file",          required_argument, 0, 0 },  /*  3 */
  {"tree_show",          no_argument,       0, 0 },  /*  4 */
  {"output_file",        required_argument, 0, 0 },  /*  5 */
  {"ptp_multi",          no_argument,       0, 0 },  /*  6 */
  {"ptp_single",         no_argument,       0, 0 },  /*  7 */
  {"outgroup",           required_argument, 0, 0 },  /*  8 */
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
  opt_quiet = 0;
  opt_ptpmulti = 0;
  opt_ptpsingle = 0;

  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) == 0)
  {
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
        opt_ptpmulti = 1;
        break;

      case 7:
        opt_ptpsingle = 1;
        break;

      case 8:
        opt_outgroup = optarg;
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
  if (opt_ptpmulti)
    commands++;
  if (opt_ptpsingle)
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
          "  --help                         display help information.\n"
          "  --version                      display version information.\n"
          "  --tree_show                    display an ASCII version of the tree.\n"
          "  --ptp_multi                    PTP style with one lambda per coalescent.\n"
          "  --ptp_single                   PTP style with single lambda for all coalescent.\n"
          "  --outgroup TAXON               In case the input tree is unrooted, use TAXON as the outgroup (default: taxon with longest branch).\n"
          "  --quiet                        only output warnings and fatal errors to stderr.\n"
          "Input and output options:\n"
          "  --tree_file FILENAME           tree file in newick format.\n"
          "  --output_file FILENAME         output file name.\n"
         );
}

void cmd_ptpmulti(bool multiple_lambda)
{
  FILE * out;

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

    rtree = utree_convert_rtree(utree, tip_count);
    utree_destroy(utree);
  }
  else
  {
    if (!opt_quiet)
      fprintf(stdout, "Loaded rooted tree...\n");
  }
  
  /* TODO: Sarah's heuristic function should be called here */
  ptp_multi_heuristic(rtree, multiple_lambda);

  if (opt_treeshow)
    rtree_show_ascii(rtree);

  if (!opt_quiet)
    fprintf(stdout, "Writing tree file...\n");

  /* export tree structure to newick string */
  char * newick = rtree_export_newick(rtree);

  /* Write newick to file */
  out = fopen(opt_outfile, "w");
  if (!out)
    fatal("Cannot write to file %s", opt_outfile);

  fprintf(out, "%s;", newick);
  fclose(out);
  free(newick);

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

  cmdline = (char *)xmalloc(len + argc + 1);
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

  if (opt_help)
  {
    cmd_help();
  }
  else if (opt_ptpmulti)
  {
    cmd_ptpmulti(true);
  }
  else if (opt_ptpsingle)
  {
    cmd_ptpmulti(false);
  }

  free(cmdline);
  return (0);
}
