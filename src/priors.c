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

// function taken from http://stackoverflow.com/questions/24294192/computing-the-binomial-coefficient-in-c
long long combi(int n, int k)
{
  long long ans=1;
  k=k>n-k?n-k:k;
  int j=1;
  for(;j<=k;j++,n--)
  {
    if(n%j==0)
    {
      ans*=n/j;
    }
    else if (ans%j==0)
    {
      ans=ans/j*n;
    }
    else
    {
      ans=(ans*n)/j;
    }
  }
  return ans;
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
  params_uniform* params = (params_uniform*) (info.params);

  double from = params->uniform_from;
  double to = params->uniform_to;
  assert(from < to);
  return 1.0 / (from - to);
}

double beta_logprior(int num_species, prior_inf info)
{
  params_beta* params = (params_beta*) (info.params);

  assert(num_species > 0);
  double alpha = params->beta_alpha;
  double beta = params->beta_beta;
  assert(alpha > 0);
  assert(beta > 0);
  double x = (double) num_species;
  return log(gamma(alpha+beta)) - log(gamma(alpha)) - log(gamma(beta))
    + log(pow(x, alpha-1)) + log(pow(1-x, beta-1));
}

double dirichlet_logprior(int num_species, prior_inf info)
{
  params_dirichlet* params = (params_dirichlet*) (info.params);

  assert(num_species > 0);
  assert(0);
  // TODO: not implemented!
  return -1;
}

double gamma_logprior(int num_species, prior_inf info)
{
  params_gamma* params = (params_gamma*) (info.params);

  assert(num_species > 0);
  double beta = params->gamma_rate;
  assert(beta > 0);
  double alpha = params->gamma_shape;
  assert(alpha > 0);

  return log(pow(beta, alpha)) - log(gamma(alpha))
    + log(pow(num_species, alpha - 1))
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

  return log(combi(n, num_species))
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

  return log(combi(k + r - 1, k)) + log(pow(1.0 - p, r)) + log(pow(p, k));
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
