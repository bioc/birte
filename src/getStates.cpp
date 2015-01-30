/**
 * @file    getStates.cpp
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
 * This file implements the C wrapper function for interfacing C++ code from R.
 */

#include <RcppArmadillo.h>
#include "BayesNetwork.h"

// TODO: Warnings/Errors, wenn Kanten nicht exisitieren
extern "C" {
	#include "getStates.h"

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
			SEXP sexpaffinitiesTF, SEXP sexpPTC_miRNA, SEXP sexpaffinitiesQ) {

		int i,j,c,r;		
		int C_cnt = INTEGER(num_conditions)[0];
		int O_cnt = INTEGER(num_mRNA)[0];
		int A_cnt = INTEGER(num_miRNA)[0];
		int T_cnt = INTEGER(num_TF)[0];
		int Q_cnt = INTEGER(n_others)[0];
		int use_miRNA_expression = INTEGER(sexp_use_miRNA_expression)[0];
		int use_Q_expression = INTEGER(sexp_use_Q_expression)[0];
		int mRNADataType = INTEGER(sexp_mRNADataType)[0];
		int miRNADataType = INTEGER(sexp_miRNADataType)[0];
		int TFexprDataType = INTEGER(sexp_TFexprDataType)[0];
		int QDataType = INTEGER(sexp_QDataType)[0];
		int thin = INTEGER(sexpthin)[0];
		int model = INTEGER(sexpmodel)[0];

		// hyperparameters
		double n0 = REAL(sexpn0)[0];     
		double alpha = REAL(sexpalpha)[0];
		double beta = REAL(sexpbeta)[0];
		double alphamiR = REAL(sexpalphamiR)[0];
		double betamiR = REAL(sexpbetamiR)[0];
		double alphaTF = REAL(sexpalphaTF)[0];
		double betaTF = REAL(sexpbetaTF)[0];
		double alphaQ = REAL(sexpalphaQ)[0];
		double betaQ = REAL(sexpbetaQ)[0];
		
		double* theta_TF = (double*) R_alloc(T_cnt, sizeof(double));
    for(i = 0; i < T_cnt; i++){
      theta_TF[i] = (double)REAL(sexptheta_TF)[i];      
    }
		double* theta_miRNA = (double*) R_alloc(A_cnt, sizeof(double));
     for(i = 0; i < A_cnt; i++){
      theta_miRNA[i] = (double)REAL(sexptheta_miRNA)[i];
    }
		double* theta_Q = NULL;
		double** L = NULL;
		int** interactions=NULL;
    int K_cnt = INTEGER(sexpK_cnt)[0];    
		if(K_cnt > 0 && Q_cnt > 0){
			Rprintf("Using informative prior for regulator states\n");     
      Rprintf("K_cnt = %i\n", K_cnt);
			L = (double**) R_alloc(K_cnt, sizeof(double*));
			for(i = 0; i < K_cnt; i++){
				L[i] = (double*) R_alloc(K_cnt, sizeof(double));
				for(j = 0; j < K_cnt; j++){
					L[i][j] = (double)REAL(sexpL)[i + j*(K_cnt)];
				}
			}
			interactions = (int**) R_alloc(Q_cnt, sizeof(int*));      
			for(i = 0; i < Q_cnt; i++){
				interactions[i] = (int*) R_alloc(2, sizeof(int));
				interactions[i][0] = (int)INTEGER(sexpinteract)[i];
				interactions[i][1] = (int)INTEGER(sexpinteract)[i + Q_cnt];        
			}
		}
    else if(K_cnt == 0 && Q_cnt > 0){      
        theta_Q = (double*) R_alloc(Q_cnt, sizeof(double));
        for(i = 0; i < Q_cnt; i++){
          theta_Q[i] = (double)REAL(sexptheta_Q)[i];
    }
    }
	
		int only_switches = INTEGER(sexponly_switches)[0];
		list<int> *S_potential_swaps = NULL;
		if((A_cnt > 0) & (!only_switches)) {
			S_potential_swaps = new list<int>[A_cnt];
			for(i=0; i<A_cnt; i++) {
				int curr_size = LENGTH(VECTOR_ELT(sexpS_potential_swaps, i));
				for(j=0; j<curr_size; j++) {
					// Adding IDs of potenital swap partners, indices change in C
					S_potential_swaps[i].push_back(INTEGER(VECTOR_ELT(sexpS_potential_swaps, i))[j]-1);
				}
		
			}
			if(S_potential_swaps == NULL)
				Rprintf("S_potential_swaps: memory allocation problem!\n");
		}

		list<int> *T_potential_swaps = NULL;
		if((T_cnt > 0) & (!only_switches)) {
			T_potential_swaps = new list<int>[T_cnt];
			for(i=0; i<T_cnt; i++) {
				int curr_size = LENGTH(VECTOR_ELT(sexpT_potential_swaps, i));
				for(j=0; j<curr_size; j++) {
					// Adding IDs of potenital swap partners, indices change in C
					T_potential_swaps[i].push_back(INTEGER(VECTOR_ELT(sexpT_potential_swaps, i))[j]-1);
				}
			}
			if(T_potential_swaps == NULL)
				Rprintf("T_potential_swaps: memory allocation problem!\n");

		}

		list<int> *Q_potential_swaps = NULL;
		if((Q_cnt > 0) & (!only_switches)) {
			Q_potential_swaps = new list<int>[Q_cnt];
			for(i=0; i<Q_cnt; i++) {
				int curr_size = LENGTH(VECTOR_ELT(sexpQ_potential_swaps, i));
				for(j=0; j<curr_size; j++) {
					// Adding IDs of potenital swap partners, indices change in C
					Q_potential_swaps[i].push_back(INTEGER(VECTOR_ELT(sexpQ_potential_swaps, i))[j]-1);
				}
			}
			if(Q_potential_swaps == NULL)
				Rprintf("Q_potential_swaps: memory allocation problem!\n");

		}
    Rprintf("Q_cnt = %i\n", Q_cnt);
		
		double **affinitiesTF = (double**) R_alloc(T_cnt, sizeof(double*));
		for(i=0; i<T_cnt; i++) {
			int curr_size = LENGTH(VECTOR_ELT(sexpaffinitiesTF, i));
			affinitiesTF[i] = (double*)R_alloc(curr_size, sizeof(double));
			for(j=0; j<curr_size; j++) {
				affinitiesTF[i][j] = REAL(VECTOR_ELT(sexpaffinitiesTF, i))[j];
			}
		}
		double **PTC_miRNA = (double **)R_alloc(A_cnt, sizeof(double*));
		for(i=0; i<A_cnt; i++) {
			int curr_size = LENGTH(VECTOR_ELT(sexpPTC_miRNA, i));
			PTC_miRNA[i] = (double*)R_alloc(curr_size, sizeof(double));
			for(j=0; j<curr_size; j++) {
				PTC_miRNA[i][j] = REAL(VECTOR_ELT(sexpPTC_miRNA, i))[j];
			}
		}
		double **affinitiesQ = (double**) R_alloc(Q_cnt, sizeof(double*));
		for(i=0; i<Q_cnt; i++) {
			int curr_size = LENGTH(VECTOR_ELT(sexpaffinitiesQ, i));
			affinitiesQ[i] = (double*)R_alloc(curr_size, sizeof(double));
			for(j=0; j<curr_size; j++) {
				affinitiesQ[i][j] = REAL(VECTOR_ELT(sexpaffinitiesQ, i))[j];
			}
		}

		int** init_S = NULL;
		if(A_cnt > 0){
			init_S = (int**)R_alloc(C_cnt, sizeof(int*));
			for(c = 0; c < C_cnt; c++){
				init_S[c] = (int*) R_alloc(A_cnt, sizeof(int));
				for(i = 0; i < A_cnt; i++){
					init_S[c][i] = (int)INTEGER(sexpinit_S)[c + i*C_cnt];
					//Rprintf("S[%i][%i] = %i ", c, i, init_S[c][i]);
				}
			}
		}
		// IDs of transcription factors
		int** init_T = NULL;
		if(T_cnt > 0){
			init_T = (int**) R_alloc(C_cnt, sizeof(int*));
			for(c = 0; c < C_cnt; c++){
				init_T[c] = (int*) R_alloc(T_cnt, sizeof(int));
				for(i = 0; i < T_cnt; i++){
					init_T[c][i] = (int)INTEGER(sexpinit_T)[c + i*C_cnt];
				}
			}
		}
		int** init_Q = NULL;
		if(Q_cnt > 0){
			init_Q = (int**) R_alloc(C_cnt, sizeof(int*));
			for(c = 0; c < C_cnt; c++){
				init_Q[c] = (int*) R_alloc(Q_cnt, sizeof(int));
				for(i = 0; i < Q_cnt; i++){
					init_Q[c][i] = (int)INTEGER(sexpinit_Q)[c + i*C_cnt];
				}
			}
		}
	      // #replicates of experiment (miRNA=0, mRNA=1) e under condition c (control=0, treated=1), access rep_cnt[e][c]
	  int **rep_cnt = (int**)R_alloc(4, sizeof(int *));
	  int k;
	  for(k = 0; k < 4; k++){
		  rep_cnt[k] = (int *)R_alloc(C_cnt, sizeof(int)); // 0=miRNA, 1=mRNA, 2=TF, 3=others
	  }
	  r = 0;
	  for(int k = 0; k < 4; k++){
			for(c = 0; c < C_cnt; c++){
				rep_cnt[k][c] = INTEGER(replicates)[r];
				//Rprintf("rep_cnt[%i][%i] = %i (first index: data type, second index: condition)\n", k, c, rep_cnt[k][c]);
				r++;
			}
	  }

		int nrep;
		double ***A = NULL;
		double *A_sigma = NULL;
		double *alpha_i0 = NULL;
		double *alpha_i = NULL;
		if(use_miRNA_expression){
			A = (double ***) R_alloc(C_cnt, sizeof(double**));
			nrep = 0;
			for(c=0; c<C_cnt; c++) {
				A[c] = (double**) R_alloc(A_cnt, sizeof(double*));
				for(i=0; i<A_cnt; i++) {
					A[c][i] = (double*) R_alloc(rep_cnt[0][c], sizeof(double));
					for(r=0; r<rep_cnt[0][c]; r++) {
						A[c][i][r] = REAL(miRNA_expr)[i+((r+nrep)*A_cnt)];
					}
				}
				nrep += rep_cnt[0][c];
			}
			alpha_i0 = (double *)R_alloc(A_cnt, sizeof(double));
			alpha_i = (double *) R_alloc(A_cnt, sizeof(double));
			A_sigma = (double *) R_alloc(A_cnt, sizeof(double));
			for(i=0; i<A_cnt; i++) {
				A_sigma[i] = REAL(miRNA_sigma)[i];
				alpha_i0[i] = REAL(sexpalpha_i0)[i];
				alpha_i[i] = REAL(sexpalpha_i)[i];
			}
		}
			// Expression under condition c = {0=control,1=treated} of mRNA/miRNA i, in replicate r
		double ***O = (double ***) R_alloc(C_cnt, sizeof(double**));
		nrep = 0;
		for(c=0; c<C_cnt; c++) {
			O[c] = (double**) R_alloc(O_cnt, sizeof(double*));
			for(i=0; i<O_cnt; i++) {
				O[c][i] = (double*) R_alloc(rep_cnt[1][c], sizeof(double));
				for(r=0; r<rep_cnt[1][c]; r++) {
					O[c][i][r] = REAL(mRNA_expr)[i+((r+nrep)*O_cnt)];
				}
			}
			nrep += rep_cnt[1][c];
		}
		double *O_sigma = NULL;
		O_sigma = (double *) R_alloc(O_cnt, sizeof(double));
		for(i=0; i<O_cnt; i++) {
			O_sigma[i] = REAL(mRNA_sigma)[i];
		}

		double*** Qdat = NULL;
		double* Q_sigma = NULL;
		double *alpha_i0Q = NULL;
		double *alpha_iQ = NULL;
		if(use_Q_expression){
			Qdat = (double ***) R_alloc(C_cnt, sizeof(double**));
			nrep = 0;
			for(c=0; c<C_cnt; c++) {
				Qdat[c] = (double**) R_alloc(Q_cnt, sizeof(double*));
				for(i=0; i<O_cnt; i++) {
					Qdat[c][i] = (double*) R_alloc(rep_cnt[3][c], sizeof(double));
					for(r=0; r<rep_cnt[3][c]; r++) {
						Qdat[c][i][r] = REAL(Q_expr)[i+((r+nrep)*Q_cnt)];
					}
				}
				nrep += rep_cnt[3][c];
			}

			Q_sigma = (double *) R_alloc(Q_cnt, sizeof(double));
			alpha_i0Q = (double *)R_alloc(Q_cnt, sizeof(double));
			alpha_iQ = (double *) R_alloc(Q_cnt, sizeof(double));
			for(i=0; i<Q_cnt; i++) {
				Q_sigma[i] = REAL(others_sigma)[i];
				alpha_i0Q[i] = REAL(sexpalpha_i0Q)[i];
				alpha_iQ[i] = REAL(sexpalpha_iQ)[i];
			}
		}

		int nTFexpr = INTEGER(sexpnTFexpr)[0];
		double ***Otf = NULL;
		double *alpha_i0TF = NULL;
		double *alpha_iTF = NULL;
		double *TF_sigma = NULL;
		//Rprintf("%d\n", nTFexpr);
		if(nTFexpr > 0) {
			// Expression for transcription factors
			Otf = (double ***) R_alloc(C_cnt, sizeof(double**));
			nrep = 0;
			for(c=0; c<C_cnt; c++) {
				Otf[c] = (double**) R_alloc(nTFexpr, sizeof(double*));
				for(i=0; i<nTFexpr; i++) {
					Otf[c][i] = (double*) R_alloc(rep_cnt[2][c], sizeof(double));
					for(r=0; r<rep_cnt[2][c]; r++) {
						Otf[c][i][r] = REAL(sexpTFexpr)[i+((r+nrep)*nTFexpr)];
						//Rprintf("%f ", Otf[c][i][r]);
					}
					//Rprintf("\n");
				}
				nrep += rep_cnt[2][c];
			}

			alpha_i0TF = (double *)R_alloc(nTFexpr, sizeof(double));
			alpha_iTF = (double *) R_alloc(nTFexpr, sizeof(double));
			TF_sigma = (double *) R_alloc(nTFexpr, sizeof(double));
			for(i=0; i<nTFexpr; i++) {
				alpha_i0TF[i] = REAL(sexpalpha_i0TF)[i];
				alpha_iTF[i] = REAL(sexpalpha_iTF)[i];
				TF_sigma[i] = REAL(sexpTF_sigma)[i];
				//Rprintf("%f %f %f\n", alpha_i0TF[i], alpha_iTF[i], TF_sigma[i]);
			}
		}

	
		list<int> *S2O = NULL;
		list<int> *SparentsOfO = NULL;
		if(A_cnt > 0) {
			// Intialization of edges between miRNAs and genes
			// Note that indices in R start at 1. Thus index-1 here
			S2O = new list<int>[A_cnt];
			SparentsOfO = new list<int>[O_cnt];
			for(i=0; i<A_cnt; i++) {
				for(j=0; j<LENGTH(VECTOR_ELT(mirTargets, i)); j++) {
					int currMirTarget = INTEGER(VECTOR_ELT(mirTargets, i))[j]-1;
					SparentsOfO[currMirTarget].push_back(i);
					S2O[i].push_back(currMirTarget);
				}
			}
			if(S2O == NULL | SparentsOfO == NULL)
				Rprintf("S2O: memory allcation problem!\n");
		}

		list<int> *T2O = NULL;
		list<int> *TparentsOfO = NULL;
		if(T_cnt > 0) {
			// Intialization of edges between TFs and genes
			T2O = new list<int>[T_cnt];
			TparentsOfO = new list<int>[O_cnt];
			for(i=0; i<T_cnt; i++) {
				for(j=0; j<LENGTH(VECTOR_ELT(TFtargets, i)); j++) {
					int currTFtarget = INTEGER(VECTOR_ELT(TFtargets, i))[j]-1;
					T2O[i].push_back(currTFtarget);
					TparentsOfO[currTFtarget].push_back(i);
				}
			}
			if(T2O == NULL | TparentsOfO == NULL)
				Rprintf("T2O: memory allcation problem!\n");

		}

		list<int> *Q2O = NULL;
		list<int> *QparentsOfO = NULL;
		if(Q_cnt > 0) {
			// Intialization of edges between additional factors and genes
			Q2O = new list<int>[Q_cnt];
			QparentsOfO = new list<int>[O_cnt];
			for(i=0; i<Q_cnt; i++) {
				for(j=0; j<LENGTH(VECTOR_ELT(QTargets, i)); j++) {
					int currQtarget = INTEGER(VECTOR_ELT(QTargets, i))[j]-1;
					Q2O[i].push_back(currQtarget);
					QparentsOfO[currQtarget].push_back(i);
				}
			}
			if(Q2O == NULL | QparentsOfO == NULL)
				Rprintf("Q2O: memory allcation problem!\n");

		}

		int niterations = INTEGER(niter)[0];
		int burnin = INTEGER(sexpburnin)[0];		
		BayesNetwork bn(C_cnt, O_cnt, A_cnt, T_cnt, Q_cnt, rep_cnt, O, A, Otf, Qdat,
				mRNADataType, miRNADataType, TFexprDataType, QDataType, nTFexpr,
				O_sigma, A_sigma, TF_sigma, Q_sigma,
				S2O, SparentsOfO, T2O, TparentsOfO, Q2O, QparentsOfO,
				alpha, beta,
				n0, alphamiR, betamiR, alphaTF, betaTF, alphaQ, betaQ,
				alpha_i0, alpha_i, alpha_i0TF, alpha_iTF, alpha_i0Q, alpha_iQ,
				model, only_switches,
				S_potential_swaps, T_potential_swaps, Q_potential_swaps,
				theta_TF, theta_miRNA, theta_Q, L, interactions,
				init_S, init_T, init_Q, affinitiesTF, PTC_miRNA, affinitiesQ);

		Rprintf("sampling ...\n");
		double *log_lik_trace = new double[niterations + burnin + 1];
    bn.MCMC(niterations, burnin, thin, log_lik_trace);
		Rprintf("finished.\n");
		
    arma::mat** mycoef = bn.getCoef();
    arma::mat map = bn.getMAP();        
    
    // Create R output
    int R = 0;
    for(c = 0; c < C_cnt; c++){
      if(rep_cnt[1][c] > R)
        R = rep_cnt[1][c];
    }
    arma::field<arma::mat> mycoefF(C_cnt, R);
    for(c = 0; c < C_cnt; c++){
      for(r = 0; r < rep_cnt[1][c]; r++)
        mycoefF(c, r) = mycoef[c][r];
    }	
    arma::mat ltrace(log_lik_trace, niterations + burnin + 1, 1);		    

		if(A_cnt > 0){
			delete[] S2O;
			delete[] SparentsOfO;
		}
		if(T_cnt > 0){
			delete[] TparentsOfO;
			delete[] T2O;
		}
		if(Q_cnt > 0){
			delete[] QparentsOfO;
			delete[] Q2O;
		}
		if(S_potential_swaps != NULL)
			delete[] S_potential_swaps;
		if(T_potential_swaps != NULL)
			delete[] T_potential_swaps;
		if(Q_potential_swaps != NULL)
				delete[] Q_potential_swaps;
    delete[] log_lik_trace;    
        		
    SEXP sexp;
		sexp = PROTECT(Rcpp::wrap(Rcpp::List::create(
							Rcpp::Named("post", bn.getPost()),
              Rcpp::Named("map", map),
							Rcpp::Named("log_lik_trace", ltrace),
							Rcpp::Named("eff_sample_size", bn.getEffSampleSize()),
              Rcpp::Named("coef", mycoefF))));              
		UNPROTECT(1);            
		return sexp;
	}
  
  /*SEXP getPPD(SEXP sexpcoef, SEXP sexpx){            
    arma::mat x = Rcpp::as<arma::mat>(sexpx);    
    arma::mat coef = Rcpp::as<arma::mat>(sexpcoef);    
    int** rep_cnt = new int*[4];
    for(int t = 0; t < 4; t++){
      rep_cnt[t] = new int[1];
    }    
    rep_cnt[1][0] = 1;
    arma::mat** acoef = new arma::mat*[1];
    acoef[0] = new arma::mat[1];
    acoef[0][0] = coef;
    BayesNetwork bn(acoef, rep_cnt);   
    arma::mat pred = bn.postPredDistr(0, 0, x);            
    delete[] rep_cnt;
    delete[] acoef;
    SEXP sexp;
  	PROTECT(sexp = Rcpp::wrap(Rcpp::List::create(
							Rcpp::Named("pred", pred))));
		UNPROTECT(1);
		return sexp;
  }*/
  
}
