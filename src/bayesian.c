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

void ptp_bayesian(rtree_t * rtree, bool multiple_lambda, double p_value,
  bool quiet, double min_br, int prior, int hyperprior)
{
  // select prior
  PRIOR_FUNC prior_function;
  prior_inf prior_info;
  int num_taxa = rtree->leaves;
  switch(prior)
  {
    case PRIOR_UNIFORM:
      prior_function = uniform_logprior;
      params_uniform* params_uni_prior = malloc(sizeof(params_uniform));
      params_uni_prior->uniform_from = 1;
      params_uni_prior->uniform_to = num_taxa;
      prior_info.params = params_uni_prior;
      break;

    case PRIOR_NEGATIVE_BINOMIAL:
      prior_function = negative_binomial_logprior;
      params_negative_binomial* params_neg_binom =
        malloc(sizeof(params_negative_binomial));
      params_neg_binom->negative_binomial_probability = 0.1337;
      params_neg_binom->negative_binomial_failures = 42;
      prior_info.params = params_neg_binom;
      break;

    case PRIOR_BINOMIAL:
      prior_function = binomial_logprior;
      params_binomial* params_binom = malloc(sizeof(params_binomial));
      params_binom->binomial_probability = 0.1337;
      params_binom->binomial_n = num_taxa;
      prior_info.params = params_binom;
      break;

    case PRIOR_GAMMA:
      prior_function = gamma_logprior;
      params_gamma* params_g = malloc(sizeof(params_gamma));
      params_g->gamma_rate = 42;
      params_g->gamma_shape = 42;
      prior_info.params = params_g;
      break;

    case PRIOR_DIRICHLET:
      prior_function = dirichlet_logprior;
      params_dirichlet* params_d = malloc(sizeof(params_dirichlet));
      prior_info.params = params_d;
      break;

    case PRIOR_BETA:
      prior_function = beta_logprior;
      params_beta* params_b = malloc(sizeof(params_beta));
      params_b->beta_alpha = 42;
      params_b->beta_beta = 42;
      prior_info.params = params_b;
      break;

    default:
      prior_function = no_logprior;
  }

  // select hyperprior
  HYPERPRIOR_FUNC hyperprior_function;
  hyperprior_inf hyperprior_info;
  switch(hyperprior)
  {
    case HYPERPRIOR_EXPONENTIAL:
      hyperprior_function = exponential_hyperprior;
      params_exponential* params_exp = malloc(sizeof(params_exponential));
      params_exp->exponential_rate = 42;
      hyperprior_info.params = params_exp;
      break;

    default:
      hyperprior_function = uniform_hyperprior;
      params_uniform* params_uni = malloc(sizeof(params_uniform));
      params_uni->uniform_from = 42;
      params_uni->uniform_to = 1337;
      hyperprior_info.params = params_uni;
  }

  free(prior_info.params);
  free(hyperprior_info.params);
}
