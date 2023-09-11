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

double IncompleteGamma(double x, double alpha, double ln_gamma_alpha)
{
   /* returns the incomplete gamma ratio I(x,alpha) where x is the upper
              limit of the integration and alpha is the shape parameter.
      returns (-1) if in error
      ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
      (1) series expansion,     if (alpha>x || x<=1)
      (2) continued fraction,   otherwise
      RATNEST FORTRAN by
      Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
      19: 285-287 (AS32)
   */
   int i;
   double p = alpha, g = ln_gamma_alpha;
   double accurate = 1e-10, overflow = 1e60;
   double factor, gin = 0, rn = 0, a = 0, b = 0, an = 0, dif = 0, term = 0, pn[6];

   if (x == 0) return (0);
   if (x < 0 || p <= 0) return (-1);

   factor = exp(p*log(x) - x - g);
   if (x > 1 && x >= p) goto l30;
   /* (1) series expansion */
   gin = 1;  term = 1;  rn = p;
l20:
   rn++;
   term *= x / rn;   gin += term;
   if (term > accurate) goto l20;
   gin *= factor / p;
   goto l50;
l30:
   /* (2) continued fraction */
   a = 1 - p;   b = a + x + 1;  term = 0;
   pn[0] = 1;  pn[1] = x;  pn[2] = x + 1;  pn[3] = x*b;
   gin = pn[2] / pn[3];
l32:
   a++;
   b += 2;
   term++;
   an = a*term;
   for (i = 0; i < 2; i++)
      pn[i + 4] = b*pn[i + 2] - an*pn[i];
   if (pn[5] == 0) goto l35;
   rn = pn[4] / pn[5];
   dif = fabs(gin - rn);
   if (dif > accurate) goto l34;
   if (dif <= accurate*rn) goto l42;
l34:
   gin = rn;
l35:
   for (i = 0; i < 4; i++) pn[i] = pn[i + 2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i = 0; i < 4; i++) pn[i] /= overflow;
   goto l32;
l42:
   gin = 1 - factor*gin;

l50:
   return (gin);
}

#define CDFGamma(x,alpha,beta) IncompleteGamma((beta)*(x),alpha,lgamma(alpha))
#define CDFChi2(x,v) CDFGamma(x,(v)/2.0,0.5)



double loglikelihood(long edge_count, double edgelen_sum)
{
  assert(edge_count >= 0);

  if (edge_count == 0 || edgelen_sum < __DBL_MIN__) return 0;

  return edge_count * (log(edge_count) - 1 - log(edgelen_sum));
}

int lrt(double nullmodel_logl, double ptp_logl, unsigned int df, double * pvalue)
{
  double diff = 2*(ptp_logl - nullmodel_logl);
  *pvalue = 1 - CDFChi2(diff,df);

  if ((*pvalue) > opt_pvalue)
    return 0;

  return 1;
}

double aic(double logl, long k, long n)
{
  if (k > 1) k++;

  return -2*logl + 2*k + (double)(2*k*(k + 1)) / (double)(n-k-1);
}
