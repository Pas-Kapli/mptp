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
  bool quiet, double min_br, int prior, int hyperprior_1, int hyperprior_2,
  int runs)
{
  assert(runs > 0);

  // select prior and hyperpriors
  PRIOR_FUNC prior_function;
  prior_inf prior_info;
  int num_taxa = rtree->leaves;

  HYPERPRIOR_FUNC hyperprior_function_1;
  hyperprior_inf hyperprior_info_1;
  HYPERPRIOR_FUNC hyperprior_function_2;
  hyperprior_inf hyperprior_info_2;
  // +/- how much in the proposals?
  double window_hyperprior_1 = 1;
  double window_hyperprior_2 = 1;

  assert(window_hyperprior_1 > 0);
  assert(window_hyperprior_2 > 0);

  int num_hyperpriors = 0;

  switch(prior)
  {
    case PRIOR_UNIFORM:
      prior_function = uniform_logprior;
      params_uniform* params_uni_prior = malloc(sizeof(params_uniform));
      params_uni_prior->uniform_from = 1;
      params_uni_prior->uniform_to = num_taxa;
      prior_info.params = params_uni_prior;
      num_hyperpriors = 0;
      assert(hyperprior_1 == HYPERPRIOR_NONE);
      assert(hyperprior_2 == HYPERPRIOR_NONE);
      break;

    case PRIOR_NEGATIVE_BINOMIAL:
      prior_function = negative_binomial_logprior;
      params_negative_binomial* params_neg_binom =
        malloc(sizeof(params_negative_binomial));
      params_neg_binom->negative_binomial_probability = 0.1337;
      params_neg_binom->negative_binomial_failures = 42;
      prior_info.params = params_neg_binom;
      num_hyperpriors = 2;
      switch(hyperprior_1)
      {
        case HYPERPRIOR_EXPONENTIAL:
          hyperprior_function_1 = exponential_hyperprior;
          params_exponential* params_exp = malloc(sizeof(params_exponential));
          params_exp->exponential_rate = 42;
          hyperprior_info_1.params = params_exp;
          break;

        case HYPERPRIOR_UNIFORM:
          hyperprior_function_1 = uniform_hyperprior;
          params_uniform* params_uni = malloc(sizeof(params_uniform));
          params_uni->uniform_from = 0;
          params_uni->uniform_to = 1;
          hyperprior_info_1.params = params_uni;
          break;

        default:
          fatal("No hyperprior for the negative binomial probability selected.\n");
      }
      switch(hyperprior_2)
      {
        case HYPERPRIOR_EXPONENTIAL:
          hyperprior_function_1 = exponential_hyperprior;
          params_exponential* params_exp = malloc(sizeof(params_exponential));
          params_exp->exponential_rate = 42;
          hyperprior_info_1.params = params_exp;
          break;

        case HYPERPRIOR_UNIFORM:
          hyperprior_function_2 = uniform_hyperprior;
          params_uniform* params_uni = malloc(sizeof(params_uniform));
          params_uni->uniform_from = 42;
          params_uni->uniform_to = 1337;
          hyperprior_info_2.params = params_uni;
          break;

        default:
          fatal("No hyperprior for the negative binomial failures selected.\n");
      }

      break;

    case PRIOR_BINOMIAL:
      prior_function = binomial_logprior;
      params_binomial* params_binom = malloc(sizeof(params_binomial));
      params_binom->binomial_probability = 0.1337;
      params_binom->binomial_n = num_taxa;
      prior_info.params = params_binom;
      num_hyperpriors = 1;
      assert(hyperprior_2 == HYPERPRIOR_NONE);
      switch(hyperprior_1)
      {
        case HYPERPRIOR_EXPONENTIAL:
          hyperprior_function_1 = exponential_hyperprior;
          params_exponential* params_exp = malloc(sizeof(params_exponential));
          params_exp->exponential_rate = 42;
          hyperprior_info_1.params = params_exp;
          break;

        case HYPERPRIOR_UNIFORM:
          hyperprior_function_1 = uniform_hyperprior;
          params_uniform* params_uni = malloc(sizeof(params_uniform));
          params_uni->uniform_from = 0;
          params_uni->uniform_to = 1;
          hyperprior_info_1.params = params_uni;
          break;

        default:
          fatal("No hyperprior for the binomial probability selected.\n");
      }
      break;

    case PRIOR_GAMMA:
      prior_function = gamma_logprior;
      params_gamma* params_g = malloc(sizeof(params_gamma));
      params_g->gamma_rate = 42;
      params_g->gamma_shape = 42;
      prior_info.params = params_g;
      num_hyperpriors = 2;
      switch(hyperprior_1)
      {
        case HYPERPRIOR_EXPONENTIAL:
          hyperprior_function_1 = exponential_hyperprior;
          params_exponential* params_exp = malloc(sizeof(params_exponential));
          params_exp->exponential_rate = 42;
          hyperprior_info_1.params = params_exp;
          break;

        case HYPERPRIOR_UNIFORM:
          hyperprior_function_1 = uniform_hyperprior;
          params_uniform* params_uni = malloc(sizeof(params_uniform));
          params_uni->uniform_from = 42;
          params_uni->uniform_to = 1337;
          hyperprior_info_1.params = params_uni;
          break;

        default:
          fatal("No hyperprior for the gamma rate selected.\n");
      }
      switch(hyperprior_2)
      {
        case HYPERPRIOR_EXPONENTIAL:
          hyperprior_function_1 = exponential_hyperprior;
          params_exponential* params_exp = malloc(sizeof(params_exponential));
          params_exp->exponential_rate = 42;
          hyperprior_info_1.params = params_exp;
          break;

        case HYPERPRIOR_UNIFORM:
          hyperprior_function_2 = uniform_hyperprior;
          params_uniform* params_uni = malloc(sizeof(params_uniform));
          params_uni->uniform_from = 42;
          params_uni->uniform_to = 1337;
          hyperprior_info_2.params = params_uni;
          break;

        default:
          fatal("No hyperprior for the gamma shape selected.\n");
      }
      break;

    case PRIOR_DIRICHLET:
      prior_function = dirichlet_logprior;
      params_dirichlet* params_d = malloc(sizeof(params_dirichlet));
      prior_info.params = params_d;
      num_hyperpriors = 0;
      assert(hyperprior_1 == HYPERPRIOR_NONE);
      assert(hyperprior_2 == HYPERPRIOR_NONE);
      break;

    case PRIOR_BETA:
      prior_function = beta_logprior;
      params_beta* params_b = malloc(sizeof(params_beta));
      params_b->beta_alpha = 42;
      params_b->beta_beta = 42;
      prior_info.params = params_b;
      num_hyperpriors = 2;
      switch(hyperprior_1)
      {
        case HYPERPRIOR_EXPONENTIAL:
          hyperprior_function_1 = exponential_hyperprior;
          params_exponential* params_exp = malloc(sizeof(params_exponential));
          params_exp->exponential_rate = 42;
          hyperprior_info_1.params = params_exp;
          break;

        case HYPERPRIOR_UNIFORM:
          hyperprior_function_1 = uniform_hyperprior;
          params_uniform* params_uni = malloc(sizeof(params_uniform));
          params_uni->uniform_from = 42;
          params_uni->uniform_to = 1337;
          hyperprior_info_1.params = params_uni;
          break;

        default:
          fatal("No hyperprior for the beta alpha selected.\n");
      }
      switch(hyperprior_2)
      {
        case HYPERPRIOR_EXPONENTIAL:
          hyperprior_function_1 = exponential_hyperprior;
          params_exponential* params_exp = malloc(sizeof(params_exponential));
          params_exp->exponential_rate = 42;
          hyperprior_info_1.params = params_exp;
          break;

        case HYPERPRIOR_UNIFORM:
          hyperprior_function_2 = uniform_hyperprior;
          params_uniform* params_uni = malloc(sizeof(params_uniform));
          params_uni->uniform_from = 42;
          params_uni->uniform_to = 1337;
          hyperprior_info_2.params = params_uni;
          break;

        default:
          fatal("No hyperprior for the beta beta selected.\n");
      }
      break;

    default:
      prior_function = no_logprior;
      assert(hyperprior_1 == HYPERPRIOR_NONE);
      assert(hyperprior_2 == HYPERPRIOR_NONE);
  }

  delimit_stats* best_solution = ptp_multi_heuristic(rtree, multiple_lambda,
    p_value, true, min_br, prior_function, prior_info);

  assert(num_hyperpriors >= 0);
  assert(num_hyperpriors <= 2);

  if (num_hyperpriors > 0)
  {
    delimit_stats* previous_solution = best_solution;
    int i;
    for (i = 0; i < runs; ++i)
    {
      // do some bayesian stuff, MCMC for the hyperpriors
      // TODO: grab old parameter value, do the proposal,
      //   decide to take the jump or not
      if (num_hyperpriors >= 1)
      {
        switch(prior)
        {
          case PRIOR_NEGATIVE_BINOMIAL: // negative binomial probability
            assert(window_hyperprior_1 <= 1);
            double old_negative_binomial_probability = ((params_negative_binomial*) (prior_info.params))->negative_binomial_probability;
            break;
          case PRIOR_BINOMIAL: // binomial probability
            assert(window_hyperprior_1 <= 1);
            double old_binomial_probability = ((params_binomial*) (prior_info.params))->binomial_probability;
            break;
          case PRIOR_GAMMA: // gamma rate
            ;
            double old_gamma_rate = ((params_gamma*) (prior_info.params))->gamma_rate;
            break;
          case PRIOR_BETA: // beta alpha
            ;
            double old_beta_alpha = ((params_beta*) (prior_info.params))->beta_alpha;
            break;
          default:
            fatal("Something is wrong with the priors (1).\n");
        }
        if (num_hyperpriors == 2)
        {
          switch(prior)
          {
            case PRIOR_NEGATIVE_BINOMIAL: // negative binomial failures
              ;
              int old_negative_binomial_failures = ((params_negative_binomial*) (prior_info.params))->negative_binomial_failures;
              break;
            case PRIOR_GAMMA: // gamma shape
              ;
              double old_gamma_shape = ((params_gamma*) (prior_info.params))->gamma_shape;
              break;
            case PRIOR_BETA: // beta beta
              ;
              double old_beta_beta = ((params_beta*) (prior_info.params))->beta_beta;
              break;
            default:
              fatal("Something is wrong with the priors (2).\n");
          }
        }
      }
    }
  }

  if (prior_info.params)
  {
    free(prior_info.params);
  }

  if (num_hyperpriors >= 1)
  {
    free(hyperprior_info_1.params);
  }

  if (num_hyperpriors == 2)
  {
    free(hyperprior_info_2.params);
  }

  free(best_solution);
}
