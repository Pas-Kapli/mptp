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

double uniform_logprior(int num_species, int num_taxa)
{
  assert(num_species > 0);
  assert(num_taxa > 0);
  return log(((double) num_species) / ((double) num_taxa));
}

double no_logprior(int num_species, int num_taxa)
{
  assert(num_species > 0);
  assert(num_taxa > 0);
  return 0;
}

/*
Nice idea, but did not work :-(

PRIOR_FUNC create_uniform_logprior(int num_taxa)
{
  assert(num_taxa > 0);
  double logprior(int num_species) {
        return uniform_logprior(num_species, num_taxa);
  }
  return logprior;
}

PRIOR_FUNC create_no_logprior(int num_taxa)
{
  assert(num_taxa > 0);
  double logprior(int num_species) {
        return no_logprior(num_species, num_taxa);
  }
  return logprior;
}
*/
