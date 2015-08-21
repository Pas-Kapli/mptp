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

double exponential_hyperprior(double x, hyperprior_inf info)
{
  params_exponential* params = (params_exponential*) (info.params);

  assert(x >= 0);
  double rate = params->exponential_rate;
  return rate * exp(-rate * x);
}

double uniform_hyperprior(double x, hyperprior_inf info)
{
  params_uniform * params = (params_uniform*) (info.params);

  double from = params->uniform_from;
  double to = params->uniform_to;
  assert(from < to);
  return 1.0 / (from - to);
}

double beta_logprior(int num_species, prior_inf info)
{
  beta_params_t * params = (beta_params_t *) (info.params);

  assert(num_species > 0);
  double alpha = params->alpha;
  double beta = params->beta;
  assert(alpha > 0);
  assert(beta > 0);
  double x = (double) num_species;
  return log(gamma(alpha+beta)) - log(gamma(alpha)) - log(gamma(beta))
    + log(pow(x, alpha-1)) + log(pow(1-x, beta-1));
}

double dirichlet_logprior(int num_species, prior_inf info)
{
  params_dirichlet * params = (params_dirichlet*) (info.params);

  assert(num_species > 0);
  assert(0);
  // TODO: not implemented!
  return -1;
}

double gamma_logpdf(int num_species, prior_inf info)
{
  gamma_params_t * params = (gamma_params_t *) (info.params);
  double beta = params->beta;
  double alpha = params->alpha;
  assert(alpha > 0 && beta > 0 && num_species > 0);

  return alpha*log(beta) - log(gamma(alpha))
    + (alpha-1)*log(num_species)
    - beta * num_species;
}

double binomial_logprior(int num_species, prior_inf info)
{
  params_binomial* params = (params_binomial*) (info.params);

  assert(num_species > 0);
  double p = params->binomial_probability;
  assert(p > 0);
  assert(p < 1);
  int n = params->binomial_n;
  assert(n > 0);

  return logcomb(n, num_species)
    + log(pow(1.0 - p, num_species))
    + log(pow(p, n - num_species));
}

double negative_binomial_logprior(int num_species, prior_inf info)
{
  params_negative_binomial* params = (params_negative_binomial*) (info.params);

  assert(num_species > 0);
  double p = params->negative_binomial_probability;
  assert(p > 0);
  assert(p < 1);
  int r = params->negative_binomial_failures;
  assert(r > 0);

  // r is the number of failures until the experiment is stopped
  int k = num_species; // k is the number of successes

  return logcomb(k + r - 1, k) + log(pow(1.0 - p, r)) + log(pow(p, k));
}

double uniform_logprior(int num_species, prior_inf info)
{
  params_uniform* params = (params_uniform*) (info.params);

  assert(num_species > 0);
  double from = params->uniform_from;
  double to = params->uniform_to;
  assert(from < to);
  return 1.0 / (from - to);
}

double no_logprior(int num_species, prior_inf info)
{
  assert(num_species > 0);
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
