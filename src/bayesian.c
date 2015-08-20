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
  PRIOR_FUNC prior_function;
  switch(prior)
  {
    case PRIOR_UNIFORM:
      prior_function = uniform_logprior;
      break;

    case PRIOR_NEGATIVE_BINOMIAL:
      prior_function = negative_binomial_logprior;
      break;

    case PRIOR_BINOMIAL:
      prior_function = binomial_logprior;
      break;

    case PRIOR_GAMMA:
      prior_function = gamma_logprior;
      break;

    case PRIOR_DIRICHLET:
      prior_function = dirichlet_logprior;
      break;

    case PRIOR_BETA:
      prior_function = beta_logprior;
      break;

    default:
      prior_function = no_logprior;
  }

  HYPERPRIOR_FUNC hyperprior_function;
  switch(hyperprior)
  {
    case HYPERPRIOR_EXPONENIAL:
      hyperprior_function = exponential_hyperprior;
      break;

    default:
      hyperprior_function = uniform_hyperprior;
  }
}
