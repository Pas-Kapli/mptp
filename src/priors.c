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

static double * logn_table;

/* precompute all log(i!) for 0 <= i <= n and store in logn_table */
void init_logn_table(unsigned int n)
{
  unsigned int i;

  logn_table = (double *)xmalloc((n+1)*sizeof(double));

  logn_table[0] = 0;
  for (i = 1; i <= n; ++i)
    logn_table[i] = logn_table[i-1] + log(i);
}

/* calculate log comb(n,k) from precomputed table */
static double logcomb(unsigned int n, unsigned int k)
{
  return (logn_table[n] - (logn_table[n-k] + logn_table[k]));
}


double prior_score(unsigned int species_count, prior_t * prior)
{
  if (!prior) return 0;

  int dist = prior->dist;

  if (dist == PRIOR_NONE)
    return 0;
  else if (dist == PRIOR_UNI)
    return uni_logpdf(species_count, (uni_params_t *)(prior->params));
  else if (dist == PRIOR_BINOMIAL)
    return bin_logpmf(species_count, (bin_params_t *)(prior->params));
  else if (dist == PRIOR_NBIN)
    return nbin_logpmf(species_count, (nbin_params_t *)(prior->params));
  else if (dist == PRIOR_BETA)
    return beta_logpdf(species_count, (beta_params_t *)(prior->params));
  else if (dist == PRIOR_GAMMA)
    return gamma_logpdf(species_count, (gamma_params_t *)(prior->params));
  else if (dist == PRIOR_EXP)
    return exp_logpdf(species_count, (exp_params_t *)(prior->params));

  /* TODO: Handle the dirichlet distribution */
    
  return 0;
}


double beta_logpdf(double x, beta_params_t * params)
{
  assert(x > 0);

  double alpha = params->alpha;
  double beta = params->beta;
  assert(params->alpha > 0);
  assert(params->beta > 0);
  return log(gamma(alpha+beta)) - log(gamma(alpha)) - log(gamma(beta))
    + log(pow(x, alpha-1)) + log(pow(1-x, beta-1));
}

double gamma_logpdf(double x, gamma_params_t * params)
{
  double beta = params->beta;
  double alpha = params->alpha;
  assert(alpha > 0 && beta > 0 && x > 0);

  return alpha*log(beta) - 
         lgamma(alpha) +
         (alpha-1)*log(x) -
         beta*x;
}


double dir_logpdf(double * x, dir_params_t * params)
{
  assert(0);
  // TODO: not implemented!
  return -1;
}

double bin_logpmf(unsigned int k, bin_params_t * params)
{
  assert(k >= 0);
  double p = params->prob;
  assert(p >= 0 && p <= 1);
  int n = params->trials;
  assert(n >= 0);

  return logcomb(n, k) +
         log(pow(1.0 - p, k)) +
         log(pow(p, n - k));
}

double nbin_logpmf(unsigned int k, nbin_params_t * params)
{
  assert(k >= 0);
  double p = params->prob;
  assert(p >= 0 && p <= 1);
  int r = params->failures;
  assert(r > 0);

  return logcomb(k + r - 1, k) + 
         log(pow(1.0 - p, r)) + 
         log(pow(p, k));
}

double uni_logpdf(double x, uni_params_t * params)
{
  assert(x >= params->min && x <= params->max);
  assert(params->min < params->max);

  return 1.0 / (params->max - params->min);
}

double exp_logpdf(double x, exp_params_t * params)
{
  assert(params->rate > 0);
  assert(x >= 0);

  if (x < 0) return -__DBL_MAX__ / 2;

  return log(params->rate) - x*(params->rate);
}
