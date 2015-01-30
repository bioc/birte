/**
 * @file    getStates.h
 * @version 0.99.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * This file contains the function definition of the wrapper function for interfacing C++ code from R.
 */

#ifndef GETSTATES_HEADER
#define GETSTATES_HEADER

#include <R_ext/Rdynload.h>

extern "C"{

/**
  * Wrapper functions for interfacing C++ code from R
  * 
  * @return posterior state inference.
  */
// [[register]]
  SEXP getStates(SEXP num_conditions, SEXP num_mRNA,SEXP num_miRNA, SEXP num_TF, SEXP n_others, SEXP replicates,
			SEXP mRNA_expr, SEXP miRNA_expr, SEXP sexpTFexpr, SEXP Q_expr, SEXP sexpnTFexpr,
			SEXP mRNA_sigma, SEXP miRNA_sigma, SEXP sexpTF_sigma, SEXP others_sigma,
			SEXP sexp_mRNADataType, SEXP sexp_miRNADataType, SEXP sexp_TFexprDataType, SEXP sexp_QDataType,
			SEXP sexp_use_miRNA_expression, SEXP sexp_use_Q_expression,
			SEXP mirTargets, SEXP TFtargets, SEXP QTargets,
			SEXP sexpalpha, SEXP sexpbeta,
			SEXP sexpn0, SEXP sexpalphamiR, SEXP sexpbetamiR, SEXP sexpalphaTF, SEXP sexpbetaTF, SEXP sexpalphaQ, SEXP sexpbetaQ,
			SEXP sexpalpha_i0, SEXP sexpalpha_i, SEXP sexpalpha_i0TF, SEXP sexpalpha_iTF, SEXP sexpalpha_i0Q, SEXP sexpalpha_iQ,
			SEXP sexpmodel, SEXP niter, SEXP sexpburnin, SEXP sexpthin, SEXP sexponly_switches,
			SEXP sexpT_potential_swaps, SEXP sexpS_potential_swaps, SEXP sexpQ_potential_swaps,
			SEXP sexptheta_TF, SEXP sexptheta_miRNA, SEXP sexptheta_Q, SEXP sexpL, SEXP sexpK_cnt, SEXP sexpinteract,
			SEXP sexpinit_S, SEXP sexpinit_T, SEXP sexpinit_Q,
			SEXP sexpaffinitiesTF, SEXP sexpPTC_miRNA, SEXP sexpaffinitiesQ);
      
  /**
   * @return posterior predictive distribution
   */
  //SEXP getPPD(SEXP sexpcoef, SEXP sexpx);
    
R_CallMethodDef callMethods[]  = {
  {"getStates", (DL_FUNC) &getStates, 59},
  {NULL, NULL, 0}
};

void
R_init_birte(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
}
	 
#endif
