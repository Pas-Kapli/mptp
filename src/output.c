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

FILE * open_file_ext(const char * extension, long seed)
{
  char * filename = NULL;
  if (opt_mcmc)
  {
    if (asprintf(&filename, "%s.%ld.%s", opt_outfile, seed, extension) == -1)
      fatal("Unable to allocate enough memory.");
  }
  else
  {
    if (asprintf(&filename, "%s.%s", opt_outfile, extension) == -1)
      fatal("Unable to allocate enough memory.");
  }

  FILE * out = xopen(filename,"w");

  free(filename);

  return out;
}

void output_info(FILE * out,
                 long method,
                 double nullmodel_logl,
                 double logl,
                 double pvalue,
                 int lrt_result,
                 rtree_t * root,
                 unsigned int species_count)
{
  fprintf(out, "Command: %s\n", cmdline);
  fprintf(out,
          "Number of edges greater than minimum branch length: %d / %d\n",
           root->edge_count,
           2 * root->leaves - 2);
  fprintf(out, "Null-model score: %.6f\n", nullmodel_logl);
  fprintf(out,
          "Best score for %s coalescent rate: %.6f\n",
          (method == PTP_METHOD_SINGLE) ?
            "single" : "multi",
          logl);
  fprintf(out, "LRT computed p-value: %.6f\n", pvalue);
  fprintf(out, "LRT: %s\n", lrt_result ? "passed" : "failed");
  fprintf(out, "Number of delimited species: %d\n", species_count);
}

void output_minbr(double minbr)
{
  FILE * fp_out;
  char * filename = NULL;

  if (asprintf(&filename, "%s.%s", opt_outfile, "txt") == -1)
    fatal("Unable to allocate enough memory.");

  fp_out = xopen(filename,"w");
  free(filename);

  fprintf(fp_out, "Command: %s\n", cmdline);
  if (minbr == 0)
    fprintf(fp_out, "Minimum branch length (--minbr) should be set to 0\n");
  else
    fprintf(fp_out, "Minimum branch length (--minbr) should be set to %f\n",minbr);
}
