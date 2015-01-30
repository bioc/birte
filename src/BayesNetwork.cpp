/**
 * @file    BayesNetwork.cpp
 * @author  Holger Froehlich
 * @version 0.99.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the free Software Foundation; either version 2 of
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
 * This file implements all functions and classes defined in the corresponding HEADER file.
 */



#include "BayesNetwork.h"

using namespace std;


BayesNetwork::BayesNetwork(){
}

BayesNetwork::BayesNetwork(arma::mat** coef, int** rep_cnt){
  this->coef = coef;  
  this->rep_cnt = rep_cnt;
}
		
BayesNetwork::BayesNetwork(int C_cnt, int O_cnt, int A_cnt, int T_cnt, int Q_cnt, int **rep_cnt,
		double ***mRNA_expression, double ***miRNA_expression, double ***Otf, double*** Q_expression,
		int mRNADataType, int miRNADataType,  int TFexprDataType, int QDataType,  int nTFexpr,
		double *mRNA_sigma, double *miRNA_sigma, double *TF_sigma, double* Q_sigma,
		list<int> *S2O, list<int> *SparentsOfO, list<int> *T2O, list<int> *TparentsOfO, list<int>* Q2O, list<int>* QparentsOfO,
		double alpha, double beta,
		double n0, double alphamiR, double betamiR, double alphaTF, double betaTF, double alphaQ, double betaQ,
		double *alpha_i0, double* alpha_i, double *alpha_i0TF, double *alpha_iTF, double* alpha_i0Q, double* alpha_iQ,
		int model, int only_switches,
		list<int> *S_potential_swaps, list<int> *T_potential_swaps, list<int>* Q_potential_swaps,
		double* theta_TF, double* theta_miRNA, double* theta_Q, double** K, int** interactions,
		int** init_S, int** init_T, int** init_Q, double** affinitiesTF, double** affinitiesmiRNA, double** affinitiesQ) {

	this->MODEL = NO_PLUG_IN;
	this->O_cnt = O_cnt;
	this->A_cnt = A_cnt;
	this->T_cnt = T_cnt;
	this->Q_cnt = Q_cnt;
	this->C_cnt = C_cnt;
	this->rep_cnt = rep_cnt;
  
  this->affinitiesTF = affinitiesTF;
  this->affinitiesmiRNA = affinitiesmiRNA;
	this->affinitiesQ = affinitiesQ;

	//this->mRNAs = mRNAs;
	//this->miRNAs = miRNAs;
	//this->TFs = TFs;

	this->A = miRNA_expression;
	this->O = mRNA_expression;
	this->Otf = Otf;
	this->Qdat = Q_expression;
	this->nTFexpr = nTFexpr;
	this->O_sigma = mRNA_sigma;
	this->A_sigma = miRNA_sigma;
	this->Q_sigma = Q_sigma;
	this->TF_sigma = TF_sigma;
	this->mRNADataType = mRNADataType;
	this->miRNADataType = miRNADataType;
	this->TFexprDataType = TFexprDataType;
	this->QDataType = QDataType;
	
	if(C_cnt > 2 && (A != NULL || Otf != NULL || Qdat != NULL)){
		Rprintf("Warning: MiRNA, TF and other data ignored. MiRNA, TF and other data can only be integrated with C_cnt = 2 conditions or relative expression changes!\n");
		A = NULL;
		Otf = NULL;
		Qdat = NULL;
	}

	this->S2O = S2O;
	this->T2O = T2O;
	this->Q2O = Q2O;
	this->SparentsOfO = SparentsOfO;
	this->TparentsOfO = TparentsOfO;
	this->QparentsOfO = QparentsOfO;

	int i,j,c,r;
	// hyperparameters
	this->alpha = alpha;
	this->beta = beta;
	this->n0 = n0;
	this->alphamiR = alphamiR;
	this->betamiR = betamiR;
	this->alphaTF = alphaTF;
	this->betaTF = betaTF;
	
	this->theta_TF = theta_TF;
	this->theta_miRNA = theta_miRNA;
	this->theta_Q = theta_Q;	
	this->K = K;	
	this->interactions = interactions;
	this->eps = new double[C_cnt];
  for(c = 0; c < C_cnt; c++)
	  eps[c] = 1;	
    
	this->alpha_i0 = alpha_i0;
	this->alpha_i = alpha_i;	
	this->alpha_i0TF = alpha_i0TF;
	this->alpha_iTF = alpha_iTF;
	this->alpha_i0Q = alpha_i0Q;
	this->alpha_iQ = alpha_iQ;
	//this->omega_miRNA = omega_miRNA; 
	//this->omega_TF = omega_TF;

  this->nselected = new int*[C_cnt];
  for(c = 0; c < C_cnt; c++){
		this->nselected[c] = new int[3];   
    this->nselected[c][0] = 0;
    this->nselected[c][1] = 0;
    this->nselected[c][2] = 0;
	}  
	// Initialization of S	
	//this->posterior_miRNA = (double **)CALLOC(C_cnt, double*);
	this->S = new int*[C_cnt];
	for(c=0; c<C_cnt; c++) {
	  this->S[c] = new int[A_cnt];
	  //this->posterior_miRNA[c] = (double *) CALLOC(A_cnt, double);
	  for(i=0; i<A_cnt; i++) {
	        S[c][i] = init_S[c][i];          
          this->nselected[c][1] += S[c][i];          
		//this->posterior_miRNA[c][i] = 0;
	  }
		//Rprintf("\n\n");
	}
	// Initialization of T
	//this->posterior_TF = (double **)CALLOC(C_cnt, double*);
	this->T =new int*[C_cnt];
	for(c=0; c<C_cnt; c++) {
	  this->T[c] = new int[T_cnt];
	 // this->posterior_TF[c] = (double *) CALLOC(T_cnt, double);
	  for(i=0; i<T_cnt; i++) {
	        T[c][i] = init_T[c][i];
           this->nselected[c][0] += T[c][i];
		//this->posterior_TF[c][i] = 0;
	  }
	}
	// Initialization of Q
	//this->posterior_Q = (double **)CALLOC(C_cnt, double*);
	this->Q = new int*[C_cnt];
	for(c=0; c<C_cnt; c++) {
	  this->Q[c] = new int[Q_cnt];
	 // this->posterior_Q[c] = (double *) CALLOC(Q_cnt, double);
	  for(i=0; i<Q_cnt; i++) {
			Q[c][i] = init_Q[c][i];
       this->nselected[c][2] += Q[c][i];      
		//this->posterior_Q[c][i] = 0;
	  }
	}

	// These lists are empty, because no TF is active
	this->T_potential_swaps = T_potential_swaps;
	T_possible_swaps = new set<int>*[C_cnt];
	T_swap_idx = new set<int>[C_cnt];
	this->S_potential_swaps = S_potential_swaps;
	S_possible_swaps = new set<int>*[C_cnt];
	S_swap_idx = new set<int>[C_cnt];
	this->Q_potential_swaps = Q_potential_swaps;
	Q_possible_swaps = new set<int>*[C_cnt];
	Q_swap_idx = new set<int>[C_cnt];
	for(c=0; c<C_cnt; c++) {		
		T_possible_swaps[c] = new set<int>[T_cnt];
		S_possible_swaps[c] = new set<int>[A_cnt];
		Q_possible_swaps[c] = new set<int>[Q_cnt];
	}
		
	this->only_switches = only_switches;

	X = new arma::mat[C_cnt];
	R = new arma::mat[C_cnt];
	Rinv = new arma::mat[C_cnt];
	U = new arma::mat[C_cnt];
	Y = new arma::colvec*[C_cnt];
	Y_var = new double*[C_cnt];
	logdet = new double[C_cnt];
	bn_term = new double[C_cnt];
	
	//targets = new list<int>[C_cnt];
	regulators = new list<pair<int,int> >[C_cnt];	
	//omega_miRNA = (double**) CALLOC(C_cnt, double*);
	//omega_TF = (double**) CALLOC(C_cnt, double*);
	//omega_Q = (double**) CALLOC(C_cnt, double*);
  coef = new arma::mat*[C_cnt];
	for(c = 0; c < C_cnt; c++){
		//omega_TF[c] = (double*) CALLOC(T_cnt, double);
		//omega_miRNA[c] = (double*) CALLOC(A_cnt, double);
		//omega_Q[c] = (double*) CALLOC(Q_cnt, double);
		Y[c] = new arma::colvec[rep_cnt[1][c]];
		Y_var[c] = new double[rep_cnt[1][c]];
    logdet[c] = 0.0;
    bn_term[c] = 0.0;
    coef[c] = new arma::mat[rep_cnt[1][c]];      
	}

	post = arma::mat(A_cnt + T_cnt + Q_cnt, C_cnt);
	post.fill(0);
	map = arma::mat(A_cnt + T_cnt + Q_cnt, C_cnt);
	map.fill(0);
}

BayesNetwork::~BayesNetwork() {
	// CFree or delete statements for all allocated objects	
	int c,i,j;
	if(A_cnt > 0){
		//CFree(miR_higher_in_condition);
		for(c = 0; c < C_cnt; c++){
			//CFree(omega_miRNA[c]);
			//CFree(posterior_miRNA[c]);
			delete[] S[c];		
		}		
		//for(i = 0; i < eff_sample_size; i++){
		//	for(c=0; c < C_cnt; c++)
		//		CFree(posterior_omega_miRNA[i][c]);
		//	CFree(posterior_omega_miRNA[i]);
		//}
		//CFree(posterior_omega_miRNA);
		//CFree(omega_miRNA);
		//CFree(posterior_miRNA);		
		if(S_possible_swaps != NULL){		
			for(c = 0; c < C_cnt; c++){
				delete [] S_possible_swaps[c];
			}
			delete[] S_possible_swaps;
			delete[] S_swap_idx;
		}
    delete[] S;
	}	
	if(T_cnt > 0){
		for(c=0; c<C_cnt; c++) {
		//	CFree(omega_TF[c]);
			//CFree(posterior_TF[c]);
			delete[] T[c];
		}
		//for(i = 0; i < eff_sample_size; i++){
		//	for(c=0; c < C_cnt; c++)
		//		CFree(posterior_omega_TF[i][c]);
		//	CFree(posterior_omega_TF[i]);
		//}
		//CFree(posterior_omega_TF);
		//CFree(omega_TF);
		//CFree(posterior_TF);
		delete[] T;
		if(T_possible_swaps != NULL) {
			for(c=0; c<C_cnt; c++) {
				delete [] T_possible_swaps[c];
			}
			delete[] T_possible_swaps;			
			delete[] T_swap_idx;
		}
	}
	if(Q_cnt > 0){
		for(c=0; c<C_cnt; c++) {
		//	CFree(omega_Q[c]);
			//CFree(posterior_Q[c]);
			delete[] Q[c];
		}
		//for(i = 0; i < eff_sample_size; i++){
		//	for(c=0; c < C_cnt; c++)
		//		CFree(posterior_omega_Q[i][c]);
		//	CFree(posterior_omega_Q[i]);
		//}
		//CFree(posterior_omega_Q);
		//CFree(omega_Q);
		//CFree(posterior_Q);
		delete[] Q;
		if(Q_possible_swaps != NULL) {
			for(c=0; c<C_cnt; c++) {
				delete [] Q_possible_swaps[c];
			}
			delete[] Q_possible_swaps;
			delete[] Q_swap_idx;      
		}
	}  
	delete[] X;
	delete[] R;
	delete[] Rinv;
	delete[] U;
	delete[] regulators;  
	delete[] eps;
	delete[] logdet;
	delete[] bn_term;
	for(c = 0; c < C_cnt; c++){
		delete[] Y[c];
		delete[] Y_var[c];
		delete[] nselected[c];
    delete[] coef[c];
	}	
	delete[] nselected;
  delete[] coef;
}

/*double BayesNetwork::logNB(double x, double mu, double phi){
	double size = 1/phi; // phi ist der SCV (dispersion) parameter aus DESeq. Dieser ist 1/size.
	return(dnbinom_mu(x, size, mu, 1));
}*/


double BayesNetwork::get_mu0(double alpha_base, double alpha_add, int c, int newstate, int basestate) {
	// for C_cnt = 2: two-state ANOVA model
	if(C_cnt == 2){
		if(c == 0) // we assume the first condition to be the reference!!!!
			return(alpha_base);
		else{ // expected log FC
			return(alpha_base + alpha_add*abs(newstate - basestate));
		}
	}
	else
		return(alpha_add*newstate);
} 

// prior for nu:=1/lambda: logarithm of exponential distribution with parameter theta
/*double BayesNetwork::logHyperPrior(double lambda, double theta){
	return(log(theta) - theta / lambda);// '/' ist RICHTIG!!!
}*/

// construct design matrix from scratch for a given target and regulator list
void BayesNetwork::DesignMatrix(int c) {
	X[c] = arma::ones<arma::mat>(O_cnt, 1); // first column contains intercept	
	pair<int,int> p, p2;	
	for(list<pair<int,int> >::iterator it = regulators[c].begin(); it != regulators[c].end(); ){
		p = *it;
		if(p.second == 1){ // miRNA
			expandDesignMatrix(p.first, p.second, c, S2O, affinitiesmiRNA, 0);			
		}
		else if(p.second == 2){ // other
			expandDesignMatrix(p.first, p.second, c, Q2O, affinitiesQ, 0);			
		}
		else if(p.second == 0){ // TF
			expandDesignMatrix(p.first, p.second, c, T2O, affinitiesTF, 0);		
		}
		else
			Rprintf("Error: p.second = %i\n", p.second);		
		it++;
	}
}

// expand design matrix X[c]
int BayesNetwork::expandDesignMatrix(int switchid, int mir, int c, list<int>* edges, double** affinities, int updateRegulators){
	int found, omega_index = 0, i, j;
	list<pair<double,int> > newcol;
	list<int>::iterator it1, it2;
	double influence;
	X[c] = arma::join_rows(X[c], arma::zeros<arma::vec>(O_cnt));
	for(it1 = edges[switchid].begin(); it1 != edges[switchid].end(); it1++){
		j = *it1; // mRNA index
		influence = affinities[switchid][omega_index];		
		X[c](j, X[c].n_cols - 1) = influence;
		omega_index++;
	}

	if(updateRegulators){
		regulators[c].push_back(make_pair(switchid, mir));
		//Rprintf("regulator (%i,%i) added\n", switchid, mir);
	}
	return(-1); // nothing was deleted
}

int BayesNetwork::shrinkDesignMatrix(int c, int switchid, int mir){
	list<int>::iterator it1;
	list<pair<int,int> >::iterator it2;
	pair<int,int> p;
	int j=1, k, del=-1;	 // j startet bei 1, weil die 1. Spalte der Intercept ist
	for(it2 = regulators[c].begin(); it2 != regulators[c].end(); it2++){
		p = *it2;
		if(p.first == switchid & p.second == mir){
			////Rprintf("Regulator (%i, %i) deleted\n", p.first, p.second);
			it2 = regulators[c].erase(it2);
			del = j;
			break;
		}
		j++;
	}
	X[c].shed_col(del);
	return(del);		
}

// update design matrix by deleting one column or adding one column
int BayesNetwork::updateDesignMatrix(int c, int switchid, int mir){
	if(mir == 1){
		if(S[c][switchid] == 1){ // S[c][switchid] was switched to 0
			return(shrinkDesignMatrix(c, switchid, mir));
		}
		else{ 			
			return(expandDesignMatrix(switchid, mir, c, S2O, affinitiesmiRNA, 1));
		}
	}
	else if(mir == 2){
		if(Q[c][switchid] == 1){ // Q[c][switchid] was switched to 0
			return(shrinkDesignMatrix(c, switchid, mir));
		}
		else{
			return(expandDesignMatrix(switchid, mir, c, Q2O, affinitiesQ, 1));
		}
	}
	else{
		if(T[c][switchid] == 1){ // T[c][switchid] was switched to 0
			return(shrinkDesignMatrix(c, switchid, mir));
		}
		else{						
			return(expandDesignMatrix(switchid, mir, c, T2O, affinitiesTF, 1));
		}
	}
}

/*
// Brand's SVD update algorithm: too slow
int BayesNetwork::svd_rank1update(int c, int add, arma::colvec& v){
	arma::mat U, V;
	if(add < 0)
		v = -v;
	arma::colvec m = sv[c].U.t() * v; // m = t(sv$u)%*%v
	arma::vec p = v - sv[c].U * m; // p = sv$u%*%m
	double Ra = sqrt(arma::dot(p, p)); // Ra = sqrt(crossprod(p))
	if(Ra < 1e-17) // no success
		return FALSE;
	arma::vec P = p / Ra;
	if(add > 0){
		arma::colvec Q = arma::zeros<arma::rowvec>(X[c].n_cols + 1); // note: q = b!
		Q(Q.n_rows - 1) = 1;

		arma::rowvec r = arma::zeros<arma::rowvec>(sv[c].sing_val.n_elem + 1);
		r(r.n_rows - 1) = Ra;
		arma::mat K = arma::join_vert(arma::join_horiz(arma::diagmat(sv[c].sing_val), m), r); // Vorteil: K ist viel kleiner als X!
		arma::svd_econ(U, sv[c].sing_val, V, K);
		sv[v].U = arma::join_horiz(sv[c].U, P) * U;
		sv[c].V = arma::join_horiz(arma::join_vert(sv[c].V, arma::zeros<arma::rowvec>(sv[c].V.n_cols)), Q) * V;
	}
	else{
		arma::colvec n = sv[c].V.row(sv[c].V.n_rows - 1).t(); //
		arma::vec q = -sv[c].V*n;
		q(V.n_rows - 1) = 1 + q(V.n_rows - 1);
		Rb = sqrt(arma::dot(q, q));
		if(Rb < 1e-17) // no success
				return FALSE;
		q = q / Rb;
		arma::mat S = arma::diagmat(sv[c].sing_val);
		arma::mat K = arma::join_vert(arma::join_horiz(S, arma::zeros<arma::mat>(S.n_rows, 1)), arma::zeros<arma::rowvec>(S.n_cols + 1));
		arma::mat T = arma::join_vert(m, Ra) * join_horiz(n, Rb);
		K = K + T;
		arma::svd_econ(U, sv[c].sing_val, V, K);
		sv[c].V = arma::join_horiz(sv[c].V.shed_row(sv[c].V.n_rows - 1), q.shed_row(q.n_elem - 1)) * V.shed_col(V.n_cols -1);
		sv[c].U = arma::join_horiz(sv[c].U, P) * U.shed_col(U.n_cols -1);
		for(int i = 0; i < sv[c].V.n_cols; i++){
			if(arma::norm(sv[c].V.col(i), 2) < 1e-10){
				sv[c].V.shed_col(i);
				sv[c].U.shed_col(i);
				sv[c].sing_val.shed_row(i);
				break;
			}
		}
		sv[c].sing_val.shed_row(sv[c].n_rows - 1);
	}
	return TRUE;
}*/

//Cholesky factor rank 1 update / downdate
/*
arma::mat& BayesNetwork::cholupdate(arma::mat& A, arma::rowvec& x, int sign){
	double r, c, s, w;
	int p = A.n_rows;
	for(int k = 0; k < p; k++){
		w = A(k,k)*A(k,k) + sign*x(k)*x(k);
		r = sqrt(w);
		c = r / A(k,k);
		s = x(k) / A(k,k);
		A(k,k) = r;
		if(k + 1 < p){
			A(k, arma::span(k+1, p-1)) = (A(k, arma::span(k+1, p-1)) + sign*s*x.subvec(k+1, p-1)) / c;
			x.subvec(k+1, p-1) = c*x.subvec(k+1, p-1) - s*A(k, arma::span(k+1, p-1));
		}
	}
	return(A);
}

// re-triangulize R[c] by Housholder reflections
arma::mat& BayesNetwork::rotate(arma::mat& A){
	arma::colvec z, uk, vk;
	double tmp;
	int m = A.n_cols;
	for(int k = 0; k < m; k++){
		z = A(arma::span(k, m-1), k);
		uk = z;
		tmp = arma::norm(z, 2);
		if(z(0) > 0)
			uk(0) += tmp;
		if(z(0) < 0)
			uk(0) -= tmp;
		uk = uk / arma::norm(uk, 2);
		vk = arma::zeros<arma::colvec>(m);
		vk(arma::span(k, m-1)) = uk;
		A = A - ((2*vk) * (vk.t() * A));
	}
	return A;
}

double BayesNetwork::updateEps(int c, double deltaeps){
	double logdet_old = logdet[c];
	double bn_term_old = bn_term[c];
	int n = R[c].n_rows;
	if(deltaeps < 0){
		R[c] = arma::chol(X[c].t() * X[c] + eps[c]*arma::eye(X[c].n_cols, X[c].n_cols));
	}
	else{
		deltaeps = sqrt(deltaeps);
		arma::rowvec v;
		for(int i = 0; i < n; i++){
			v = arma::zeros<arma::rowvec>(n);
			v(i) = deltaeps;
			R[c] = cholupdate(R[c], v, 1);
		}
	}
	R[c] = trimatu(R[c]);
	Rinv[c] = arma::trimatu(arma::solve(R[c], arma::eye(R[c].n_rows, R[c].n_cols)));
	logdet[c] = 2*arma::sum(log(abs(R[c].diag())));
	U[c] = X[c]*Rinv[c];
	arma::colvec v2;
	updateEvidence2(c, FALSE, 0, v2);
	return(0.5*rep_cnt[1][c]*(logdet_old - logdet[c]) + (bn_term_old - bn_term[c]));
}*/

//update design matrix, targets and regulators sets and precompute terms that change in marginal log-likelihood
void BayesNetwork::updateEvidence(int c, int switchid, int mir, bool full_update){
	int cholupdate = FALSE;
	int add = 0, ncol, del;
	arma::colvec v;
	arma::mat Xcov;	
	if(regulators[c].size() > 0 && !full_update){
		ncol = X[c].n_cols;
		del = updateDesignMatrix(c, switchid, mir);
		//Rprintf("updated design matrix! X[%i]: %i x %i, del=%i\n", c, X[c].n_rows, X[c].n_cols, del);
		cholupdate = TRUE;
		add = X[c].n_cols - ncol;
		if(add > 0)
			v = X[c].col(X[c].n_cols - 1); // added column
		else
			add = -del;
	}
	else{
		if(switchid != -1)
			regulators[c].push_back(make_pair(switchid, mir));
		DesignMatrix(c);
		//Rprintf("de novo construction! X[%i]: %i x %i, switchid=%i, regulators[%i].size = %i\n", c, X[c].n_rows, X[c].n_cols, switchid, c, regulators[c].size());    
		Xcov = X[c].t() * X[c];		
    bool suc = false;
    int i = 0;
    while(!suc && i < 10){
      Xcov = Xcov + eps[c]*arma::eye(Xcov.n_rows, Xcov.n_cols);
	    suc = arma::chol(R[c] , Xcov); // it's the *upper* triangular matrix    
      i++;
    }     
    if(!suc){
      Xcov.print("Xcov possibly not PD!");
    }
		R[c] = trimatu(R[c]);
		Rinv[c] = arma::trimatu(arma::solve(R[c], arma::eye(R[c].n_rows, R[c].n_cols)));
		logdet[c] = 2*arma::sum(log(R[c].diag()));
		U[c] = X[c]*Rinv[c];
	}
	updateEvidence2(c, cholupdate, add, v);
}

void BayesNetwork::updateEvidence2(int c, int chol_update, int add, arma::colvec& v){
	arma::rowvec t, bt;
	// Cholesky factors: R, R^(-1) and U updates
	if(chol_update){
		if(add > 0){
			//Rprintf("column added! (#add = %i)", add);
			arma::rowvec vt = v.t();
			bt = vt * U[c];
			arma::vec b = bt.t();
			arma::mat deltamat =  arma::sqrt(vt*v - bt*b + eps[c]);
			double delta = arma::as_scalar(deltamat);
			arma::mat Z = arma::zeros<arma::rowvec>(R[c].n_cols);
			R[c] = arma::trimatu(arma::join_vert(arma::join_horiz(R[c], b), arma::join_horiz(Z, deltamat)));
			Rinv[c] = arma::trimatu(arma::join_vert(arma::join_horiz(Rinv[c], -Rinv[c]*b / delta), arma::join_horiz(Z, 1/deltamat)));
			logdet[c] += 2*log(fabs(R[c](R[c].n_rows - 1, R[c].n_cols - 1)));
		}
		else{
			add = -add;
			//Rprintf("column %i, deleted: Rinv[%i]: %i x %i; X[%i]: %i x %i   ", add, c, Rinv[c].n_rows, Rinv[c].n_cols, c, X[c].n_rows, X[c].n_cols);
			if(add != R[c].n_cols){
				/*R[c] = join_horiz(R[c], R[c].col(add)); // shift deleted column to the end
				R[c].shed_col(add);
				R[c] = rotate(R[c]); // re-triangulize
				R[c].shed_col(R[c].n_cols - 1);
				R[c].shed_row(R[c].n_rows - 1);*/
				arma::mat Xcov = R[c].t() * R[c];
				Xcov.shed_col(add);
				Xcov.shed_row(add);
				R[c] = arma::chol(Xcov); // it's the *upper* triangular matrix
				R[c] = arma::trimatu(R[c]);
				Rinv[c] = arma::trimatu(arma::solve(R[c], arma::eye(R[c].n_rows, R[c].n_cols))); //nur O(n^2) ==> Rank 1 update lohnt sich nicht
				logdet[c] = 2*arma::sum(log(arma::abs(R[c].diag())));
			}
			else{
        logdet[c] -= 2*log(fabs(R[c](R[c].n_rows - 1, R[c].n_cols - 1)));
				R[c].shed_col(R[c].n_cols - 1);
				R[c].shed_row(R[c].n_rows - 1);
				Rinv[c].shed_col(Rinv[c].n_cols - 1);
				Rinv[c].shed_row(Rinv[c].n_rows - 1);        
			}
		}
		//Rprintf("-->Rinv[%i]: %i x %i; X[%i]: %i x %i\n", c, Rinv[c].n_rows, Rinv[c].n_cols, c, X[c].n_rows, X[c].n_cols);
		U[c] = X[c]*Rinv[c]; // same dimension as X
	}
	bn_term[c] = 0;
	for(int r = 0; r < rep_cnt[1][c]; r++){
		t = Y[c][r].t() * U[c];
		bn_term[c] += log(beta + 0.5*(Y_var[c][r] - arma::as_scalar(t * t.t())));
	}
	bn_term[c] *= (alpha + 0.5*O_cnt);
}		


void BayesNetwork::update_swaps(list<int>* potential_swaps, set<int>** possible_swaps, set<int>* swap_idx, int** states, int switchid, int old_state, int condition) {
	// Remove all entries of switchid in states[c][regulator_k] with regulator_k in possible_swaps[condition][switchid]
	int curr;
	for(set<int>::iterator it = possible_swaps[condition][switchid].begin(); it != possible_swaps[condition][switchid].end(); it++) {
		curr = *it;
		possible_swaps[condition][curr].erase(switchid);
		if(possible_swaps[condition][curr].size() == 0)
			swap_idx[condition].erase(curr);
	}
	
	// Determine all new possible swap partners j and add them to possible_swaps[condition][switchid]. Also add switchid to possible_swaps[condition][j] for all possible partners j.
	possible_swaps[condition][switchid].clear();
	int added = 0, j;
	for(list<int>::iterator it = potential_swaps[switchid].begin(); it != potential_swaps[switchid].end(); it++) {
		j = *it;
		// add new entry j in possible_swaps[condition][switchid] and new entry switchid in possible_swaps[condition][j]
		if(states[condition][j] == old_state) {
			possible_swaps[condition][switchid].insert(j);
			if(!added){
				swap_idx[condition].insert(switchid);
				added = 1;
			}
			possible_swaps[condition][j].insert(switchid);
			swap_idx[condition].insert(j);
		}
	}
	if(possible_swaps[condition][switchid].size() == 0)
		swap_idx[condition].erase(switchid);
}

double BayesNetwork::swap_states(int** states, int swapid1, int mir, set<int>* swap_idx, list<int>* potential_swaps, set<int>** possible_swaps, double old_likelihood){
	// condition is determined
	int pos_swaps = 0, condition, condition1, condition2, counter;
	double delta_statePrior;
	set<int>::iterator it;
	for(condition = 0; condition < C_cnt; condition++){
		if(swapid1 >= swap_idx[condition].size()) {
				swapid1 = swapid1 - swap_idx[condition].size();
		}
		else
			break;
	}
  if(condition == C_cnt)
    condition = C_cnt - 1;
	// get element #swapid1  
	counter = 0;
	for(it = swap_idx[condition].begin(); it != swap_idx[condition].end(); it++){
		if(counter == swapid1)
			break;
		counter++;
	}
	swapid1 = *it;
	int r = getrand(possible_swaps[condition][swapid1].size()); // get random swap partner  
	// get element #r
	counter = 0;
	for(it =  possible_swaps[condition][swapid1].begin(); it != possible_swaps[condition][swapid1].end(); it++){
		if(counter == r)
			break;
		counter++;
	}
	int swapid2 = *it;
	condition1 = condition2 = condition;
  if(states[condition1][swapid1] == states[condition2][swapid2])
    return(0);

	// save current variables
	arma::mat Xold(X[condition]), R_old(R[condition]), Rinv_old(Rinv[condition]), U_old(U[condition]);
	list<pair<int,int> > regulators_old = regulators[condition];
	double logdet_old = logdet[condition];
	double bn_term_old = bn_term[condition];
	//double Ksum_old = Ksum[condition];

	// For calculation of new log-likelihood, two switches are performed
	double delta_loglik = 0;
	updateEvidence(condition1, swapid1, mir, false); // update design matrix for first switch
	states[condition1][swapid1] = -(states[condition1][swapid1]-1); // perform 1rst switch
	updateEvidence(condition2, swapid2, mir, false); // update design matrix and simulate effect on mRNA data of BOTH switches together
	delta_loglik = doSwitch(logdet_old, bn_term_old, FALSE, -(states[condition1][swapid1]-1), swapid1, condition1, mir); // simulate 1rst switch and influence on miRNA and TF data
	delta_loglik += doSwitch(logdet_old, bn_term_old, TRUE, states[condition2][swapid2], swapid2, condition2, mir);// total log-likelihood difference: Effekt auf mRNA-Daten wird eingeschlossen!

	delta_statePrior = deltaStatePrior(-(states[condition1][swapid1]-1), mir, swapid1) + deltaStatePrior(states[condition2][swapid2], mir, swapid2);	

	// get random number
	double u = 0;
	while(u == 0) {
		u = unif_rand();
	}			
	if(log(u) <= delta_loglik + delta_statePrior) { // Remark: Prior and neighbourhood do not change
		// updating possible swaps and changing states
		update_swaps(potential_swaps, possible_swaps, swap_idx, states, swapid1, -(states[condition1][swapid1]-1), condition1); // Achtung: der ALTE Zustand wird übergeben!
		update_swaps(potential_swaps, possible_swaps, swap_idx, states, swapid2, states[condition2][swapid2], condition2);
		states[condition2][swapid2] = -(states[condition2][swapid2]-1);
		//updateOmegas(condition);
	}
	else{
		// set variables back to initial state
		X[condition] = Xold;	;
		logdet[condition] = logdet_old;
		bn_term[condition] = bn_term_old;
		R[condition] = R_old;
		Rinv[condition] = Rinv_old;
		U[condition] = U_old;
		regulators[condition] = regulators_old;
		states[condition1][swapid1] = -(states[condition1][swapid1]-1); // set state of swapid1 back to initial state
		delta_loglik = 0;
		//Ksum[condition] = Ksum_old;
	}			
	return(delta_loglik);				
}

double BayesNetwork::switch_states(int** states, int switchid, int mir, set<int>* swap_idx, list<int>* potential_swaps, set<int>** possible_swaps, double old_likelihood){
	int condition;
	int n_Tpos_swaps_total = 0, n_Spos_swaps_total = 0, n_Qpos_swaps_total = 0;
	for(condition = 0; condition < C_cnt; condition++){
		n_Spos_swaps_total += S_swap_idx[condition].size();
		n_Tpos_swaps_total += T_swap_idx[condition].size();
		n_Qpos_swaps_total += Q_swap_idx[condition].size();
		if(mir == 1){ // Grund: Falls switchid nicht in der ersten Condition ist
			if(switchid >= A_cnt)
				switchid = switchid - A_cnt;
			else
				break;
		}
		else if(mir == 2){
			if(switchid >= Q_cnt)
				switchid = switchid - Q_cnt;
			else
				break;
		}
		else{
			if(switchid >= T_cnt)
				switchid = switchid - T_cnt;
			else
				break;
		}
	}	
	/*if(mir && miR_higher_in_condition != NULL){
		condition = abs(miR_higher_in_condition[switchid]);
	}
	else if(!mir && TF_higher_in_condition != NULL){
		condition = abs(TF_higher_in_condition[switchid]);
	}*/
	// save current variables
	//arma::colvec beta_old(beta_coef[condition]);
	arma::mat Xold(X[condition]), R_old(R[condition]), Rinv_old(Rinv[condition]), U_old(U[condition]);
	//list<int> targets_old = targets[condition];
	list<pair<int,int> > regulators_old = regulators[condition];
	double bn_term_old = bn_term[condition];
	double logdet_old = logdet[condition];
	//double Ksum_old = Ksum[condition];
  
	// simulate switch
	updateEvidence(condition, switchid, mir, false); // update design matrix, etc.
	double delta_loglik = doSwitch(logdet_old, bn_term_old, 1, states[condition][switchid], switchid, condition, mir);// log-likelihood difference  

	int old_neighbourhood , new_neighbourhood;
	old_neighbourhood = C_cnt*(A_cnt + T_cnt + Q_cnt);
	new_neighbourhood = old_neighbourhood;
	if(!this->only_switches) { // bei einem switch ändert sich die Anzahl der möglichen swaps
		old_neighbourhood += n_Tpos_swaps_total + n_Spos_swaps_total;
		update_swaps(potential_swaps, possible_swaps, swap_idx, states, switchid, states[condition][switchid], condition);
		n_Spos_swaps_total = 0;
		n_Tpos_swaps_total = 0;
		for(int c = 0; c < C_cnt; c++){
			n_Spos_swaps_total += S_swap_idx[c].size();
			n_Tpos_swaps_total += T_swap_idx[c].size();
			n_Qpos_swaps_total += Q_swap_idx[c].size();
		}
		new_neighbourhood = new_neighbourhood + n_Tpos_swaps_total + n_Spos_swaps_total + n_Qpos_swaps_total;
	}
	double delta_prior = deltaStatePrior(states[condition][switchid], mir, switchid); 
	// get random number
	double u = 0;
	while(u == 0) {
		u = unif_rand();
	}
	// switch is kept if:
	if(log(u) <= delta_loglik + delta_prior + log((double)new_neighbourhood/(double)old_neighbourhood)) {
		//Rprintf("regulators[%i].size=%i, logdet[%i]=%g, delta_loglik=%g, delta_prior =%g\n", condition, regulators[condition].size(), condition, logdet[condition], delta_loglik, delta_prior);
		states[condition][switchid] = -(states[condition][switchid]-1); // !!!  perform switch
		if(states[condition][switchid] == 1) // vorher 0, jetzt 1!
			nselected[condition][mir]++;
		else // vorher 1, jetzt 0!
			nselected[condition][mir]--;
    if(nselected[condition][mir] < 0)
      Rprintf("mir=%i, condition = %i, state[%i] = %i, #regulators = %i\n", mir, condition, switchid, states[condition][switchid], regulators[condition].size());
		//Rprintf("move accepted: old state = %i, new state = %i!\n\n", -(states[condition][switchid]-1), states[condition][switchid]);
	}
	else{
		if(!only_switches){
			update_swaps(potential_swaps, possible_swaps, swap_idx, states, switchid, -(states[condition][switchid]-1), condition);
		}
		X[condition] = Xold;
		logdet[condition] = logdet_old;
		bn_term[condition] = bn_term_old;
		R[condition] = R_old;
		Rinv[condition] = Rinv_old;
		U[condition] = U_old;
		regulators[condition] = regulators_old;
		delta_loglik = 0;		
	}
	return(delta_loglik);
}

void BayesNetwork::MCMC(long niter, long burnin, int thin, double* log_lik_trace) {
  	
	// stores values of log likelihood, if moves are accepted, Inf otherwise	
	long l, i;
	int j, k, c, r;	
	int mir, tf;
	list<int>::iterator start, stop;
	for(l=0; l<niter+burnin+1; l++) {
		log_lik_trace[l] = (double)INFINITY;
	}
	long eff_samples = niter / thin + 1;
	/*this->posterior_omega_TF = (double***) CALLOC(eff_samples, double***);
	this->posterior_omega_miRNA = (double***) CALLOC(eff_samples, double***);
	this->posterior_omega_Q = (double***) CALLOC(eff_samples, double***);
	for(i = 0; i < eff_samples; i++){
		this->posterior_omega_TF[i] = (double**) CALLOC(C_cnt, double**);
		this->posterior_omega_miRNA[i] = (double**) CALLOC(C_cnt, double**);
		this->posterior_omega_Q[i] = (double**) CALLOC(C_cnt, double**);
		for(c = 0; c < C_cnt; c++){
				this->posterior_omega_TF[i][c] = (double*) CALLOC(T_cnt, double);
				this->posterior_omega_miRNA[i][c] = (double*) CALLOC(A_cnt, double);
				this->posterior_omega_Q[i][c] = (double*) CALLOC(Q_cnt, double);
		}
	}*/
	// initialize regulator list and swaps
	for(c = 0; c < C_cnt; c++){
		for(j = 0; j < A_cnt; j++){
			if(S[c][j] > 0){
				regulators[c].push_back(make_pair(j, 1));
				if(!only_switches)
					update_swaps(S_potential_swaps, S_possible_swaps, S_swap_idx, S, j, 0, c);
			}
		}
		for(j = 0; j < T_cnt; j++){
			if(T[c][j] > 0){
				regulators[c].push_back(make_pair(j, 0));
				if(!only_switches)
					update_swaps(T_potential_swaps, T_possible_swaps, T_swap_idx, T, j, 0, c);
			}
		}
		for(j = 0; j < Q_cnt; j++){
			if(Q[c][j] > 0){
				regulators[c].push_back(make_pair(j, 2));
				if(!only_switches)
					update_swaps(Q_potential_swaps, Q_possible_swaps, Q_swap_idx, Q, j, 0, c);
			}
		}

		// store gene expression data in better format
		for(r = 0; r < rep_cnt[1][c]; r++){
			Y[c][r] = arma::zeros<arma::colvec>(O_cnt);
			for(j = 0; j < O_cnt; j++){
				Y[c][r](j) = O[c][j][r];
			}      
      coef[c][r] = arma::zeros<arma::mat>(A_cnt + T_cnt + Q_cnt + 1, ceil(niter / thin));
		}
		// compute design matrix and compute some initial values for marginal likelihood calculation
		for(int r = 0; r < rep_cnt[1][c]; r++){
			Y_var[c][r] = arma::as_scalar(Y[c][r].t() * Y[c][r]);
		}
		updateEvidence(c, -1, -1, true);
		Rprintf("initial X[%i]: %i x %i\n", c, X[c].n_rows, X[c].n_cols);
	}

   	double log_lik = likelihood();
   	Rprintf("initial (marginal) log-likelihood = %g\n", log_lik);
   	log_lik_trace[0] = log_lik;
	// set seed for random number generation (for MCMC-moves)
	GetRNGstate();	
	// sum of TFs and miRNAs
	int TFmiRNA_nswitch = C_cnt*(this->A_cnt  + this->T_cnt + this->Q_cnt); // number of possible switch operations (remains always the same);
	
	double delta_loglik, u, delta_prior, eps_new, eps_old, delta_eps, best_loglik=-INFINITY;	
	int nall_possible_ops, n_Tpos_swaps_total=0, n_Spos_swaps_total=0, n_Qpos_swaps_total=0, make_move, make_move_eps=0, cond=0;
  int swapid, switchid;
  arma::uvec idx;
	pair<int,int> p, p2;

	k = 0;
  for(i=0; i<niter+burnin; i++) { // MAIN LOOP
		R_CheckUserInterrupt();

		n_Tpos_swaps_total = 0;
		n_Spos_swaps_total = 0;
		n_Qpos_swaps_total = 0;
		for(c = 0; c < C_cnt; c++){
			if(T_cnt > 0)
				n_Tpos_swaps_total += T_swap_idx[c].size();
			if(A_cnt > 0)
				n_Spos_swaps_total += S_swap_idx[c].size();
			if(Q_cnt > 0)
				n_Qpos_swaps_total += Q_swap_idx[c].size();
			//Rprintf("T_swap_idx[%i].size = %i, S_swap_idx[%i].size=%i, Q_swap_idx[%i]=%i\n", c, T_swap_idx[c].size(), c, S_swap_idx[c].size(), c, Q_swap_idx[c].size());
		}
		nall_possible_ops = n_Tpos_swaps_total + n_Spos_swaps_total + n_Qpos_swaps_total;
		if(nall_possible_ops > 0 && !only_switches)
			make_move = getrand(TFmiRNA_nswitch + nall_possible_ops); // switch or swap with probability in dependency on number of possible switches and swaps
		else
			make_move = 0; // switch!  
		if(i % 100 == 0 && i > 0){ // suggest new parameters for state prior
			log_lik = likelihood();			
			Rprintf("(marginal) log-likelihood = %g, ", log_lik);			
			for(c = 0; c < C_cnt; c++){
				Rprintf("#TFs[%i] = %i, #miRNAs[%i] = %i, #others/interactions[%i] = %i ", c, nselected[c][0], c, nselected[c][1], c, nselected[c][2]);
				Rprintf("#regulators[%i] = %i ", c, regulators[c].size());
			}
			Rprintf("\n");

//			if(i % 1000 == 0){
//				make_move_eps = unif_rand() % C_cnt; // choose one eps
//				eps_new = pow(2, log2(eps[make_move_eps]) + rnorm(0, sqrt(0.01))) + 1e-7;
//				u = 0;
//				while(u) {
//						u = ((double)(unif_rand() % 100000001)) /((double)100000000);
//				}
//				delta_prior = logHyperPrior(eps_new, hyperprior_theta) - logHyperPrior(eps[make_move_eps], hyperprior_theta);
//				if(log(u) <= delta_prior){
//					eps_old = eps[make_move_eps];
//					delta_eps = eps_old - eps_new;
//					eps[make_move_eps] = eps_new;
//					Rprintf("==> new eps[%i] = %g\n", make_move_eps, eps[make_move_eps]);
//
//					// update of Cholesky factors R, Rinv and re-calculation of U, logdet; update der log-likelihood!
//					if(regulators[make_move_eps].size() > 0)
//						log_lik += updateEps(make_move_eps, delta_eps);
//				}
//			}
//			if(K != NULL && i%500 == 0){
//				cond = unif_rand() % C_cnt;
//				lambda_new = pow(2, log2(this->lambda[cond]) + rnorm(0, sqrt(0.5))) + 1e-7;
//				u = 0;
//				while(u == 0) {
//						u = ((double)(unif_rand() % 100000001)) /((double)100000000);
//				}
//				delta_prior = logHyperPrior(lambda_new, hyperprior_theta) - logHyperPrior(this->lambda[cond], hyperprior_theta) + deltaStatePrior2(this->lambda[cond], lambda_new, cond);
//        //Rprintf("hyperprior_theta = %g, lambda_old =%g, lambda_new = %g, logHyperPrior(lambda_new)=%g, logHyperPrior(lambda_old)=%g, delta_prior = %g\n", hyperprior_theta, this->lambda[cond], lambda_new, logHyperPrior(lambda_new, hyperprior_theta), logHyperPrior(this->lambda[cond], hyperprior_theta), delta_prior);
//				if(log(u) <= delta_prior){
//					lambda_old = lambda[cond];
//					this->lambda[cond] = lambda_new;
//					Rprintf("==> new lambda[%i] = %g\n", cond, this->lambda[cond]);
//				}
//			}
		}	

		if(make_move >= TFmiRNA_nswitch) { // swap
			make_move = getrand(nall_possible_ops);
			swapid = make_move;
			// determine swap move ID
			if(swapid >= n_Tpos_swaps_total + n_Spos_swaps_total && Q_cnt > 0) { // two Qs are swapped
				//Rprintf("Q swap\n");        
				swapid -= n_Tpos_swaps_total + n_Spos_swaps_total;         
				delta_loglik = swap_states(Q, swapid, 2, Q_swap_idx, Q_potential_swaps, Q_possible_swaps, log_lik);
			}
			else if(swapid >= n_Tpos_swaps_total && swapid < n_Tpos_swaps_total + n_Spos_swaps_total && A_cnt > 0) { // two miRNAs are swapped
				//Rprintf("miRNA swap\n");
				swapid -= n_Tpos_swaps_total;							
				delta_loglik = swap_states(S, swapid, 1, S_swap_idx, S_potential_swaps, S_possible_swaps, log_lik);
			} 
			else if(T_cnt > 0){ // two TFs are swapped
				//Rprintf("TF swap\n");
				delta_loglik = swap_states(T, swapid, 0, T_swap_idx, T_potential_swaps, T_possible_swaps, log_lik);
			}
			else{
				Rprintf("Error: No swap operation possible!");
			}
		} // swap
		else { // switch
			// Get random TF or miRNA for switch operation
			make_move = getrand(TFmiRNA_nswitch);
			switchid = make_move;	
			if(switchid >= C_cnt*(A_cnt + T_cnt) && Q_cnt > 0) { // Q is switched
				//Rprintf("Q switch\n");
				switchid -= C_cnt*(A_cnt + T_cnt);
				delta_loglik = switch_states(Q, switchid, 2, Q_swap_idx, Q_potential_swaps, Q_possible_swaps, log_lik);
			}
			else if(switchid >= C_cnt*A_cnt && switchid < C_cnt*(A_cnt + T_cnt) && T_cnt > 0) { // TF is switched
				//Rprintf("TF switch\n");
				switchid -= C_cnt*A_cnt;
				delta_loglik = switch_states(T, switchid, 0, T_swap_idx, T_potential_swaps, T_possible_swaps, log_lik);
			}
			// miRNA is switched
			else if(A_cnt > 0) {
				//Rprintf("miRNA switch");
				delta_loglik = switch_states(S, switchid, 1, S_swap_idx, S_potential_swaps, S_possible_swaps, log_lik);
			}
			else{
				Rprintf("Error: No switch operation possible!");
			}
		}// switch
		log_lik += delta_loglik;
		log_lik_trace[i+1] = log_lik;	    
		if(log_lik_trace[i + 1] > 1e100)
			Rprintf("Warning: log-likelihood > 1e100!\n");
		if(isnan(log_lik))
			Rprintf("Warning: log-likelihood is NA!\n");
	
		// update of samples from posterior
		if((i >= burnin) && (i%thin == 0)) {      		        
			for(c=0; c < C_cnt; c++) {    	
        if(log_lik > best_loglik)
			    map.col(c) = arma::zeros<arma::mat>(A_cnt + T_cnt + Q_cnt,1);
        idx = arma::zeros<arma::uvec>(regulators[c].size()+1); // first term is the intercept!
        j = 1;
				for(list<pair<int,int> >::iterator it = regulators[c].begin(); it != regulators[c].end(); it++) {          
					p = *it;                      
					if(p.second == 2){ // other
						p.first += T_cnt + A_cnt;					
					}
					else if(p.second == 1){ // mir
						p.first += T_cnt;						
					}					
					post(p.first, c) += 1;
					if(log_lik > best_loglik)                      
						map(p.first, c) = 1;					
          idx(j) = p.first + 1;
          j++;
				}       
        for(int r = 0; r < rep_cnt[1][c]; r++){           
            arma::mat tmp = Rinv[c]*Rinv[c].t()*X[c].t()*Y[c][r];                
            coef[c][r].submat(idx, k*arma::ones<arma::uvec>(1)) = Rinv[c]*Rinv[c].t()*X[c].t()*Y[c][r];   
        }        
			}
      if(log_lik > best_loglik)
        best_loglik = log_lik;
			k++;
		}
	} // MAIN LOOP	
	eff_sample_size = k;
	Rprintf("Effective sample size = %i\n", k);
  post = post / (double) k;  
	PutRNGstate();	
}


// compute log likelihood
double BayesNetwork::likelihood(){
	double mu0, x, log_lik = 0;
	int c, i, j, r;
	double bn; 
	double lgan = lgamma(alpha + 0.5*O_cnt);
	double lgalpha = lgamma(alpha);
	double z;
	for(c=0; c<C_cnt; c++) {
		// sum for miRNA replicates
    if(A != NULL && A_cnt > 0){
			for(i=0; i<this->A_cnt; i++) {
				mu0 = get_mu0(alpha_i0[i], alpha_i[i], c, S[c][i], S[0][i]);
				for(r=0; r<this->rep_cnt[0][c]; r++) {
					if(!isnan(A[c][i][r])) { // check for NaN
						x = (A[c][i][r] - mu0);
						if(MODEL == ALL_PLUG_IN) {
							log_lik -= x*x/A_sigma[i];
						}
						else{ // mu0 = empirical mean (=> empirical Bayes estimate!!!) viel zu konservativ ==> penalisiert jedes mu0 != 0!					
							log_lik -= (0.5 + alphamiR)*log(1+1/(2*betamiR)*x*x);
						}
					}
				}
			}
		}
		
		// sum for mRNA data ==> marginal log-likelihood
		log_lik += rep_cnt[1][c]*(0.5*O_cnt*log(eps[c]) + alpha*log(beta) + lgan - lgalpha - 0.5*logdet[c]) - bn_term[c];

		// sum for TF replicates
		if((nTFexpr > 0) && Otf != NULL) {
			for(i=0; i<this->nTFexpr; i++) {
				mu0 =  get_mu0(alpha_i0TF[i], alpha_iTF[i], c, T[c][i], T[0][i]);
				for(r=0; r<this->rep_cnt[2][c]; r++) {
					if(!isnan(Otf[c][i][r])) { // check for NaN
						x = Otf[c][i][r] - mu0;
						if(MODEL == ALL_PLUG_IN) {
							log_lik -= x*x/TF_sigma[i];
						}
						else if(MODEL == NO_PLUG_IN){						
							log_lik -= (0.5 + alphaTF)*log(1+1/(2*betaTF)*x*x);
						}
					}
				}
			}
		}

		// sum for replicates of other data (Q)
		if(Qdat != NULL && Q_cnt > 0){
			for(i=0; i<this->Q_cnt; i++) {
				mu0 = get_mu0(alpha_i0Q[i], alpha_iQ[i], c, Q[c][i], Q[0][i]);
				for(r=0; r<this->rep_cnt[3][c]; r++) {
					if(!isnan(Qdat[c][i][r])) { // check for NaN
						x = Qdat[c][i][r] - mu0;
						if(MODEL == ALL_PLUG_IN){
							log_lik -= x*x/Q_sigma[i];
						}
						else{							
							log_lik -= (0.5 + alphaQ)*log(1+1/(2*betaQ)*x*x);
						}
					}
				}
			}
		}
	}
	if(isnan(log_lik)) {
		Rprintf("Error: log_lik is NA!\n");
		return(0);
	}
	return(log_lik);
}

// compute likelihood difference before and after switch
double BayesNetwork::doSwitch(double logdet_old, double bn_term_old, int usemRNA, int oldstate, int switchid, int condition, int doMir){
	int r, i;
	double log_lik_diff = 0, mu0_before, new_mu0, x, y, bold, bnew;
	// index of regulated mRNA
	
	/*if(O_mu_old != NULL){
		for(i = 0; i < O_cnt; i++) {
			mu0_before = O_mu_old[i];
			new_mu0 = O_mu[condition][i];
			for(r=0; r<this->rep_cnt[1][condition]; r++) {
				// calculate marginal log likelihood difference for switch
				if(!isnan(O[condition][i][r])) { // check for NaN
					if(MODEL == ALL_PLUG_IN) {
						if(this->mRNADataType == ARRAY){
							x = mu0_before - O[condition][i][r];
							y = new_mu0 - O[condition][i][r];
							log_lik_diff += x*x / O_predvar_old(i) - y*y / O_predvar[condition](i);
						}
						else
							log_lik_diff += logNB(O[condition][i][r], new_mu0, O_sigma[i]) - logNB(O[condition][i][r], mu0_before, O_sigma[i]);
					}
					else if(MODEL == NO_PLUG_IN) {
						if(this->mRNADataType != ARRAY){
							Rprintf("Model %i not implemented for RNAseq data, yet!\n", MODEL);
							return(0);
						}
						log_lik_diff += (0.5 + alpha) * log( ( 1 + 1/(2*beta) * pow(O[condition][i][r]-mu0_before, 2))  / ( 1 + 1/(2*beta) * pow(O[condition][i][r]-new_mu0, 2)));
					}
				}
			}
		}
	}*/
	if(usemRNA){
		log_lik_diff += 0.5*rep_cnt[1][condition]*(logdet_old - logdet[condition]) + (bn_term_old - bn_term[condition]);
	}

	// if TF is switched, an additional term has to be evaluated
	if((doMir == 0) && (Otf != NULL) && (nTFexpr > 0) && (switchid < nTFexpr)) { // Annahme: TFs mit Daten stehen am Anfang!!!
		i = switchid; // ID of TF
		mu0_before = get_mu0(alpha_i0TF[switchid], alpha_iTF[switchid], condition, oldstate, T[0][switchid]);
		new_mu0 = get_mu0(alpha_i0TF[switchid], alpha_iTF[switchid], condition, -(oldstate-1), T[0][switchid]);
		for(r=0; r<this->rep_cnt[2][condition]; r++) {
			if(!isnan(Otf[condition][i][r])) { // check for NaN
				x = mu0_before - Otf[condition][i][r];
				y = new_mu0 - Otf[condition][i][r];        
				if(MODEL == ALL_PLUG_IN) {
					log_lik_diff += (x*x - y*y)/TF_sigma[i];
				}
				else if(MODEL == NO_PLUG_IN){
				/*	x = (Otf[condition][i][r] + mu0_before*n0) / (1 + n0);
					y = (Otf[condition][i][r] + new_mu0*n0) / (1 + n0);
					bold = 0.5*n0*mu0_before*mu0_before + 0.5*Otf[condition][i][r]*Otf[condition][i][r] + betaTF;
					bnew = 0.5*n0*new_mu0*new_mu0 + 0.5*Otf[condition][i][r]*Otf[condition][i][r] + betaTF;
					log_lik_diff +=  x*x - y*y + (alphaTF + 0.5)*(log(bold) - log(bnew));*/
					bold = -(alphaTF + 0.5) * log(1 + 1/(2*betaTF)*x*x);
					bnew = -(alphaTF + 0.5) * log(1 + 1/(2*betaTF)*y*y);          
					log_lik_diff +=  bnew - bold;          
          //Rprintf("i=%i, x = %g, y=%g, switched to = %i, alphaiTF=%g, otf=%g, new_mu0=%g, log-lik-diff = %g\n", i, x, y, -(oldstate-1), alpha_iTF[switchid], Otf[condition][i][r], new_mu0, bnew - bold);
				}        
			}
	  }
	}

	// if miRNA is switched, an additional term has to be evaluated
    if((doMir == 1) && (A != NULL)) { // für miRNA-Daten, falls nicht die vollständige likelihood gewünscht ist
		// i := ID of miRNA
		i = switchid;
		mu0_before = get_mu0(alpha_i0[switchid], alpha_i[switchid], condition, oldstate, S[0][switchid]);
		new_mu0 = get_mu0(alpha_i0[switchid], alpha_i[switchid], condition, -(oldstate-1), S[0][switchid]);
		for(r=0; r<this->rep_cnt[0][condition]; r++) {
			if(!isnan(A[condition][i][r])) { // check for NaN
				x = mu0_before - A[condition][i][r];
				y = new_mu0 - A[condition][i][r];
				if(MODEL == ALL_PLUG_IN) {
					log_lik_diff += (x*x - y*y)/A_sigma[i];
				}
				else if(MODEL == NO_PLUG_IN){
					/*x = (A[condition][i][r] + mu0_before*n0) / (1 + n0);
					y = (A[condition][i][r] + new_mu0*n0) / (1 + n0);
					bold = 0.5*n0*mu0_before*mu0_before + 0.5*A[condition][i][r]*A[condition][i][r] + betamiR;
					bnew = 0.5*n0*new_mu0*new_mu0 + 0.5*A[condition][i][r]*A[condition][i][r] + betamiR;
					log_lik_diff +=  x*x - y*y + (alphamiR + 0.5)*(log(bold) - log(bnew));*/
					bold = -(alphamiR + 0.5) * log(1 + 1/(2*betamiR)*x*x);
					bnew = -(alphamiR + 0.5) * log(1 + 1/(2*betamiR)*y*y);
					log_lik_diff +=  bnew - bold;
				}
			}
	  }
	}

    if((doMir == 2) && (Qdat != NULL)) { // für Q-Daten, falls nicht die vollständige likelihood gewünscht ist
    		// i := ID of miRNA
    		i = switchid;
    		mu0_before = get_mu0(alpha_i0Q[switchid], alpha_iQ[switchid], condition, oldstate, Q[0][switchid]);
    		new_mu0 = get_mu0(alpha_i0Q[switchid], alpha_iQ[switchid], condition, -(oldstate-1), Q[0][switchid]);
    		for(r=0; r<this->rep_cnt[0][condition]; r++) {
    			if(!isnan(Qdat[condition][i][r])) { // check for NaN
    				x = mu0_before - Qdat[condition][i][r];
    				y = new_mu0 - Qdat[condition][i][r];
    				if(MODEL == ALL_PLUG_IN) {
    					log_lik_diff += (x*x - y*y)/Q_sigma[i];
    				}
    				else if(MODEL == NO_PLUG_IN){
    					/*x = (Qdat[condition][i][r] + mu0_before*n0) / (1 + n0);
    					y = (Qdat[condition][i][r] + new_mu0*n0) / (1 + n0);
    					bold = 0.5*n0*mu0_before*mu0_before + 0.5*Qdat[condition][i][r]*Qdat[condition][i][r] + betaQ;
    					bnew = 0.5*n0*new_mu0*new_mu0 + 0.5*Qdat[condition][i][r]*Qdat[condition][i][r] + betaQ;
    					log_lik_diff +=  x*x - y*y + (alphaQ + 0.5)*(log(bold) - log(bnew));*/
    					bold = -(alphaQ + 0.5) * log(1 + 1/(2*betaQ)*x*x);
    					bnew = -(alphaQ + 0.5) * log(1 + 1/(2*betaQ)*y*y);
    					log_lik_diff += bnew - bold;
    				}
    			}
    	  }
    }   
	if(isnan(log_lik_diff)){
		Rprintf("Warning: lok_lik_diff (doSwitch) is NA!\n");
		return(0);
	}

	return log_lik_diff;
}

/*// log-prior
double BayesNetwork::statePrior(int c){
  double logprior = 0.0, theta;
  int i;
  for(i = 0; i < T_cnt; i++)
    logprior += T[c][i]*log(theta_TF[i] + 1e-10) + (1 - T[c][i]) * log(1 - theta_TF[i] + 1e-10);  
  for(i = 0; i < A_cnt; i++)
    logprior += S[c][i]*log(theta_miRNA[i] + 1e-10) + (1 - S[c][i]) * log(1 - theta_miRNA[i] + 1e-10);
  for(i = 0; i < Q_cnt; i++){    
    if(theta_Q != NULL)
      theta = theta_Q[i];
    else
      theta = K[interactions[i][0]][interactions[i][1]];
    logprior += Q[c][i]*log(theta + 1e-10) + (1 - Q[c][i]) * log(1 - theta + 1e-10);
  }      
  return logprior;
}*/

// log-prior difference
double BayesNetwork::deltaStatePrior(int oldState, int doMir, int switched_idx) { // returns DIFFERENCE of Prior between two state switches
	double deltaPrior;	
  double theta;
  if(doMir == 0)
    theta = theta_TF[switched_idx];
  else if(doMir == 1)
    theta = theta_miRNA[switched_idx];
  else{
    if(theta_Q != NULL)
      theta = theta_Q[switched_idx];
    else
      theta = K[interactions[switched_idx][0]][interactions[switched_idx][1]];
  }
  deltaPrior = oldState*log(1 - theta + 1e-10) + (1 - oldState)*log(theta + 1e-10);
	return deltaPrior;
}

/*double BayesNetwork::deltaStatePrior2(double old_lambda, double new_lambda, int cond){
	double deltaPrior = 0.0;
	if(K != NULL){
		double totalPrior = 0.0;
		for(int i = 0; i < Q_cnt; i++){
			totalPrior = totalPrior + fabs(Q[cond][i] - K[interactions[i][0]][interactions[i][1]]);     
		}    
		deltaPrior = totalPrior * (old_lambda - new_lambda);
	}  
	return deltaPrior;
}*/

/*int** BayesNetwork::getS() {
  return S;
}

int** BayesNetwork::getT() {
  return T;
}

int** BayesNetwork::getQ() {
  return Q;
}*/

/*double** BayesNetwork::getPostS() {
  return posterior_miRNA;
}		

double** BayesNetwork::getPostT() {
  return posterior_TF;
}

double** BayesNetwork::getPostQ() {
  return posterior_Q;
}*/

arma::mat& BayesNetwork::getPost(){
	return post;
}

arma::mat& BayesNetwork::getMAP(){    
  return map;
}

long BayesNetwork::getEffSampleSize(){
	return eff_sample_size;
}

/*
void BayesNetwork::copy2states(int c){
    int i;
    for(i = 0; i < T_cnt; i++)
      T[c][i] = 0;
    for(i = 0; i < A_cnt; i++)
      S[c][i] = 0;
    for(i = 0; i < Q_cnt; i++)
      Q[c][i] = 0;
    for(list<pair<int,int> >::iterator it = regulators[c].begin(); it != regulators[c].end(); it++) {
  			p = *it;
        if(p.second == 0) // TF
          T[c][p.first] = 1;
        else if(p.second == 1) // miRNA
          S[c][p.first] = 1;
        else
          Q[c][p.first] = 1;
    }
}*/

arma::mat** BayesNetwork::getCoef(){
  return this->coef;
}

/*
arma::mat& BayesNetwork::postPredDistr(int cond, int r, arma::mat& x){    
  arma::mat pred = arma::zeros<arma::mat>(x.n_rows, 2);      
  arma::mat predtmp;  
  if(cond < C_cnt && r < rep_cnt[1][cond]){
    predtmp = x * coef[cond][r];
    pred.col(1) = arma::mean(predtmp, 0);    
    pred.col(2) = arma::stddev(predtmp, 0, 0);      
  }  
  return pred;
}*/

int BayesNetwork::getrand(int n){
  int r = (int)floor(unif_rand() * n);
  if(r == n)
    r = n - 1;
  return r;
}

