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

double loglikelihood(int edge_count, double edgelen_sum)
{
  assert(edge_count >= 0);

  if (edge_count == 0 || edgelen_sum == 0) return 0;
  
  return edge_count * (log(edge_count) - 1 - log(edgelen_sum));
}

int lrt(double nullmodel_logl, double ptp_logl, unsigned int df, double * pvalue)
{
#ifdef HAVE_LIBGSL
  double diff = 2*(ptp_logl - nullmodel_logl);

  /* http://docs.scipy.org/doc/scipy/reference/generated/scipy.special.chdtr.html */
  *pvalue = 1 - gsl_cdf_chisq_P(diff,df);

  if ((*pvalue) > opt_pvalue)
    return 0;
#endif

  return 1;
}
