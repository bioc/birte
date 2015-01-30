/**
 * @file    BayesNetwork.h
 * @author  Holger Froehlich
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
 * BayesNetwork contains the miRNA-/TF-target interaction graph as an adjacency list, 
 * expression values (in a matrix) and (hyper-)parameters for the Markov-Chain-Monte-Carlo (MCMC) sampling.
 */

#ifndef BayesNetwork_HEADER
#define BayesNetwork_HEADER
	 
#define ARRAY 0
#define RNAseq 1
#define ALL_PLUG_IN 1
#define NO_PLUG_IN 2

#define CFree(x) if(x != NULL) Free(x);
#define CALLOC(nbytes, type) Calloc(nbytes, type);

#include <armadillo>
#include <new>
#include <list>
#include <utility>
#include <set>
#include <R.h>
#include <Rdefines.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>


using namespace std;

/*struct SVD{
	arma::mat U, V;
	arma::colvec sing_val;
};*/

class BayesNetwork
{
	protected:

		int MODEL; //!<@brief Two models are available. One uses exact parameters (only known, when data is simulated), the other uses conjugate priors (normal usage). See also class constructor.

		int mRNADataType; //!<@brief type of mRNA data (0 = microarray, 1 = RNAseq)
		int miRNADataType;//!<@brief type of miRNA data (0 = microarray, 1 = RNAseq)
		int TFexprDataType; //!<@brief type of TF expression data (0 = microarray, 1 = RNAseq)
		int QDataType; //!<@brief type of other expression data (0 = microarray, 1 = RNAseq)
				
		double *O_sigma; //!<@brief Variances / dispersion parameters for mRNA expression values.
		double *A_sigma; //!<@brief Variances / dispersion parameters for miRNA expression values.
		double *TF_sigma; //!<@brief Variances / dispersion parameters for TF expression values.
		double *Q_sigma; //!<@brief Variances / dispersion parameters for other expression values.
		
		// Vertices of the graph
		double ***O; //!<@brief mRNA expression values under condition c in replicate r, access: O[c][i][r]
		double ***A; //!<@brief miRNA expression values under condition c in replicate r, access: A[c][i][r]
		double ***Otf; //!<@brief TF expression values under condition c in replicate r, access: Otf[c][i][r]
		double ***Qdat; //!<@brief other expression values under condition c in replicate r, access: Qdat[c][i][r]
		int **S; //!<@brief miRNA states
		int **T; //!<@brief TF states
		int **Q; //!<@brief states for other regulators		
		int nTFexpr;
		
		
		//double **O_mu; //!<@brief PREDICTED mean expression values of mRNAs. Values are updated after each MCMC move.

		// Implementation of edges using adjacency lists
		list<int> *S2O; //!<@brief Adjacency list for edges from miRNA interaction graph. (S are miRNA vertices, O mRNA vertices).
		list<int> *T2O; //!<@brief Adjacency list for edges from TF interaction graph. (t are TF vertices, O mRNA vertices).
		list<int> *Q2O; //!<@brief Adjacency list for edges from Q interaction graph. (t are Q vertices, O mRNA vertices).
		list<int> *TparentsOfO; //!<@brief Adjacency list, which holds for a given mRNA its regulating TFs.
		list<int> *SparentsOfO; //!<@brief Adjacency list, which holds for a given mRNA its regulating miRNAs.
		list<int> *QparentsOfO; //!<@brief Adjacency list, which holds for a given mRNA its regulating miRNAs.

		int **rep_cnt; //!<@brief Number of replicates for experiment e (e=0 codes for miRNA measurements, e=1 codes for mRNA measurements, 2=TF measurements, 3=other measurements) under condition c (control=0, treated=1). Access as follows: rep_cnt[e][c].

		int O_cnt; //!<@brief number of mRNAs
		int A_cnt; //!<@brief number of miRNAs
		int T_cnt; //!<@brief number of TFs
		int Q_cnt; //!<@brief number of other regulators (Q)
		int C_cnt; //!<@brief number of conditions

		// hyperparameters
		double alpha; //!@brief parameter for conjugate Gamma prior for noise precision
		double beta; //!@brief parameter for conjugate Gamma prior for noise precision
		double n0;//!<@brief Hyperparameter n0. Needed calculation of Marginal distribution with unkown mean and unknown variance.  (normal distrubtion for miRNA expression values)
		double alphamiR;//!<@brief Hyperparameter alpha (for Gamma distributions = conjugate Prior for precision).
		double betamiR;//!<@brief Hyperparameter beta (for Gamma distributions = conjugate Prior for precision).
		double alphaTF;
		double betaTF;
		double alphaQ;
		double betaQ;
		double *alpha_i0;//!<@brief expression levels of miRNAs under reference condition. alpha_i0 is needed for the calculation for the hyperparameter mu0 (see function get_mu0).
		double* alpha_i;//!<@brief fold change for miRNAs. alpha_i is needed for the calculation for the hyperparameter mu0 (see function get_mu0).
		double *alpha_i0TF;
		double *alpha_iTF;
		double *alpha_i0Q;
		double *alpha_iQ;    		
		
		int** nselected; //!<@brief number of active regulators per regulator type (0 = TF, 1 = miRNA, 2 = Q)	
		arma::mat post; //!<@brief marginal posteriors of regulator activities.	
		arma::mat map; //!<@brief MAP regulator activities

		list<int> *T_potential_swaps;//!<@brief Contains for every TF k all potentially possible swap partners (access: T_potential_swaps[k]).
		set<int> **T_possible_swaps;//!<@brief Contains for every TF k all currently possible swap partners (access: T_potential_swaps[k]).		
		set<int>* T_swap_idx; //!<@brief index of TFs that can be swapped in each condition

		list<int> *S_potential_swaps;//!<@brief Contains for every miRNA j all potentially possible swap partners (access: S_potential_swaps[j]).
		set<int> **S_possible_swaps;//!<@brief Contains for every miRNA j all currently  possible swap partners (access: S_potential_swaps[j]).
		set<int>* S_swap_idx; //!<@brief index of miRNAs that can be swapped in each condition		

		list<int> *Q_potential_swaps;//!<@brief Contains for every Q j all potentially possible swap partners (access: Q_potential_swaps[j]).
		set<int> **Q_possible_swaps;//!<@brief Contains for every miRNA j all currently  possible swap partners (access: Q_potential_swaps[j]).
		set<int>* Q_swap_idx; //!<@brief index of Qs that can be swapped in each condition

		int only_switches;//!<@brief If only_switches equals 1, then MCMC moves only execute switches.
		

		double* theta_TF; //!<@brief probabilities to be active (TFs)
		double* theta_miRNA;//!<@brief probabilities to be active (miRNAs)
		double* theta_Q;//!<@brief probabilities to be active (Q)	
		double** K; //!<@brief probabilities of interactions (symmetric matrix)
		int** interactions; //!<@brief index of features participating in 2-way interaction terms (Q_cnt x 2 matrix)			
		double* eps; //!<@brief ridge parameter on (regularization constant for weights) for each condition		

		double** affinitiesTF; //!<@brief sequence binding affinities of TFs to target genes (access: affinitiesTF[TF][target])
		double** affinitiesmiRNA; //!<@brief sequence binding affinities of miRNA to target genes (access: affinitiesmiRNA[mIRNA][target])
		double** affinitiesQ; //!<@brief affinities of Q to target genes (access: affinitiesQ[mIRNA][target])

		arma::mat* X; // design matrix
		arma::mat* R; // Cholesky factor of Xcov (upper triangular matrix)
		arma::mat* Rinv; // inverse Cholesky factor of Xcov
		arma::mat* U; // a matrix square root of Hat (saved, because of more efficient update of R)
		double* logdet; // log determinant of regularized covariance matrix for each condition
		double* bn_term; // a_n * sum_r(log b_nr) (needed for marginal log-likelihood difference)
		double** Y_var; // variances of observation vectors for each condition and replicate
		arma::colvec** Y; // same as O, but stored in a matrix for convenience
		list<pair<int,int> >* regulators; // current list of regulators for each condition. Format: (regulator id, isMiRNA)    
    
    arma::mat** coef; // expected values of regression coefficients (#regulators x #samples) for each condition and replicate

		long eff_sample_size; // effective sample size
	public:

		/**
        	   * Constructors for BayesNetwork that sets all necessary parameters to the class' attributes.
        	   * 
        	   */
		BayesNetwork();
    
    BayesNetwork(arma::mat** coef, int** rep_cnt);

		BayesNetwork(int C_cnt, int O_cnt, int A_cnt, int T_cnt, int Q_cnt, int **rep_cnt,
				double ***mRNA_expression, double ***miRNA_expression, double ***Otf, double*** Q_expression,
				int mRNADataType, int miRNADataType,  int TFexprDataType, int QDataType,  int nTFexpr,
				double *mRNA_sigma, double *miRNA_sigma, double *TF_sigma, double* Q_sigma,
				list<int> *S2O, list<int> *SparentsOfO, list<int> *T2O, list<int> *TparentsOfO, list<int>* Q2O, list<int>* QparentsOfO,
				double alpha, double beta,
				double n0, double alphamiR, double betamiR, double alphaTF, double betaTF, double alphaQ, double betaQ,
				double *alpha_i0, double* alpha_i, double *alpha_i0TF, double *alpha_iTF, double* alpha_i0Q, double* alpha_iQ,
				int model, int only_switches,
				list<int> *S_potential_swaps, list<int> *T_potential_swaps, list<int>* Q_potential_swaps,
				double* theta_TF, double* theta_miRNA, double* theta_Q,  double** K, int** interactions,
				int** init_S, int** init_T, int** init_Q, double** affinitiesTF, double** affinitiesmiRNA, double** affinitiesQ);

		 /**
        	   * Destructor for BayesNetwork. Sets previously allocated memory of all class attributes free.
        	   *  
        	   */
		virtual ~BayesNetwork();


/** 
		  * Performs Markov-Chain-Monte-Carlo (MCMC) sampling.
		  * 
		  * @param niter #iterations after burnin
		  * @param burnin #iterations for burnin
		  * @param thin thinning of Markov chain: only use every thin's sample for posterior computation
      * @param memory of length (niter + nburnin + 1), in which to write log-likelihood values of sampled models 
		  * @return array containing the log-likelihood for all accepted steps/moves. Steps that were not accepted contain INF.
		  */
		virtual void MCMC(long niter, long burnin, int thin, double* log_lik_trace);

/**
		
		/**
		  * Simple Getter for posterior_miRNA.
		  * 
		  * @return returns posterior_miRNA.
		  */
	//	virtual double** getPostS();


		/**
		  * Simple Getter for posterior_TF.
		  * 
		  * @return returns posterior_TF.
		  */
	//	virtual double** getPostT();

		/**
		  * Simple Getter for posterior_Q.
		  *
		  * @return returns posterior_Q.
		  */
	//	virtual double** getPostQ();
  
		/**
		 * @return effective sample size
		 */
		virtual long getEffSampleSize();

	public:
		
		/**
		  * Calculates mean expression value mu_0, following the one-way ANOVA model: 
		  * 
		  * @param alpha_base base expression in condition 0
	  	  * @param alpha_add expected shift of mean expression
	  	  * @param c condition
	  	  * @param newstate hypothetical new state of X[c][i]
	  	  * @param basestate base state X[0][i]
		  * @return m_0 Mean expression value @f$ \mu_0 = \alpha_{i_0}  + \alpha_i * |s_{i_1} - s_{i_2}|  @f$
		 */
		virtual double get_mu0(double alpha_base, double alpha_add, int c, int newstate, int basestate);

		/**
		  * @return log-likelihood of data
		 */
		virtual double likelihood();

		/**
		  * Calculates difference in the log likelihood according to a switch
		  * 
		  * @param logdet_old old log-determinant
		  * @param bn_term_old old bn_term
		  * @param usemRNA should mRNA data be used
		  * @param oldstate old state of regulator
		  * @param switchid ID for miRNA or TF, which will be switched
		  * @param condition condition to be switched
		  * @param doMir Defines if miRNA (1) or TF (0) is switched
		  * @return returns difference in log-likelihood after switch.
		  */
		virtual double doSwitch(double logdet_old, double bn_term_old, int usemRNA, int oldstate, int switchid, int condition, int doMir);

		/**
		  * Function calculates difference in the log Prior between before and after switch
		  * 
		  * @param oldState state before switch.		  
		  * @param doMir indicates if miRNA is switched (1), another regulator (2) or a TF (0).
		  * @param switched_idx index of switched regulator		  
		  * @return returns log(Prior-probability)
		  */
		virtual double deltaStatePrior(int oldState, int doMir, int switched_idx);

		/**
		 * Calculate difference in log Prior after update of lambda
		 */
		//virtual double deltaStatePrior2(double old_lambda, double new_lambda, int cond);		

				
		//virtual double logHyperPrior(double lambda, double theta);
		
		/**
		expand design matrix X by one column (i.e. one regulator)
		@param switchid index of switched regulator
		@param mir miRNA (1) or not (0)
		@param c condition
		@param edges S2O or T2O per condition
		@param affinities regulator affinities for each target
		@param updateRegulators Should the regulator list be updated within the method or not?
		@return number of new columns of X
		*/
		virtual int expandDesignMatrix(int switchid, int mir, int c, list<int>* edges, double** affinities, int updateRegulators);

		/**
		remove one column of design matrix X (i.e. one regulator)
		@param c condition
		@param switchid index of switched regulator
		@param mir miRNA (1) or nor (0)
		@return index of deleted column in X
		*/
		virtual int shrinkDesignMatrix(int c, int switchid, int mir);

		/**
		update design matrix: calls expandDesignMatrix or shrinkDesignMatrix
		@param c condition
		@param switchid index of switched regulator
		@param passes return arguments of expandDesignMatrix and shrinkDesignMatrix
		*/
		virtual int updateDesignMatrix(int c, int switchid, int mir);

		/**
		compute design matrix from scratch
		@param c condition
		*/
		virtual void DesignMatrix(int c);


		/**
		 * Cholesky factor rank 1 update: R'R + sign*x*x'
		 * @param R matrix
		 * @param x
		 * @param sign +1: update, -1: downdate
		 * @return updated matrix R
		 */
		//arma::mat& cholupdate(arma::mat& R, arma::rowvec& x, int sign);

		/**
		 * Re-triangulize R[c] and Rinv[c] via Housholder reflections: equivalent to estimating R-matrix in a QR decomposition
		 * @param A matrix
		 * @return updated matrix A
		 */
		//arma::mat& rotate(arma::mat& A);

		/**
		 * after update of epsilon: update Cholesky factors R, Rinv and re-calculation of U, logdet; update of log-likelihood!
		 * @param c condition
		 * @param delta_eps: delta eps
		 * @return delta log likelihood
		 */
		//double updateEps(int c, double delta_eps);

		/**
		update model evidence after regulator switch: update design matrix, logdet and bn_term by calling updateEvidence2
		@param c condition
		@param switchid index of switched regulator
		@param mir miRNA (1) or TF (0) or something else (2)		
    @param full_update: reconstruct design matrix from scratch or not
		*/
		virtual void updateEvidence(int c, int switchid, int mir, bool full_update);

		/**
		update logdet and bn_term
		@param c condition
		@param chol_update update Cholesky factors or not?
		@param add number of added/removed columns (depending on sign)
		@param v column vector added to X
		*/
		virtual void updateEvidence2(int c, int chol_update, int add, arma::colvec& v);
		

		/** 
		  * When a regulator switches its activity, the possible swaps changes for the next move. This function updates the possible swap partners after a switch.
		  * 
		  @param potential_swaps T_potential_swaps or S_potential_swaps
		  @param possible_swaps S_possible_swaps or T_possible_swaps
		  @param swap_idx S_swap_idx or T_swap_idx
		  * @param switchid ID of TF, that has been switched.
		  * @param old_state state (0 or 1) of the TF before the switch.  
		  * @param condition condition (0 or 1), which has been switched.  
		  */
		virtual void update_swaps(list<int>* potential_swaps, set<int>** possible_swaps, set<int>* swap_idx, int** states, int switchid, int old_state, int condition);


		/**
		MCMC move: swap the states of two regulator, calls updateWeights
		@param states S or T
		@param swapid index of first picked regulator
		@param mir miRNA (1) or not (0)
		@param swap_idx S_swap_idx or T_swap_idx
	      	@param potential_swaps T_potential_swaps or S_potential_swaps
		@param possible_swaps S_possible_swaps or T_possible_swaps
		@param old_likelihood old log-likelihood#
		@param full_update rank 1 update or full update of SVD?
		@return difference in log-likelihood, if swap is accepted, 0 otherwise
		*/
		virtual double swap_states(int** states, int swapid, int mir, set<int>* swap_idx, list<int>* potential_swaps, set<int>** possible_swaps, double old_likelihood);

		/**
		MCMC move: switch the state of a regulator, calls update_swaps, updateWeights
		@param states S or T
		@param activeStates activeTFs or activeMiRNAs
		@param switchid index of switched regulator
		@param mir miRNA (1) or not (0)
		@param swap_idx S_swap_idx or T_swap_idx
		 @param potential_swaps T_potential_swaps or S_potential_swaps
		@param possible_swaps S_possible_swaps or T_possible_swaps
		@param old_likelihood old log-likelihood
		@param full_update rank 1 update or full update of SVD?
		@return difference in log-likelihood, if swap is accepted, 0 otherwise
		*/
		virtual double switch_states(int** states, int switchid, int mir, set<int>* swap_idx, list<int>* potential_swaps, set<int>** possible_swaps, double old_likelihood);
      
    
    /**
     * @return each matrix: expectation values of sampled regression coefficients (#coef x #samples after thinning). Whole matrix: #conditions x #replicates
     */
    virtual arma::mat** getCoef();
    
    /**
		 * @return marginal posterior activity probabilities (#regulators x #conditions)
		 */
		virtual arma::mat& getPost();
    
    /**
     * @return configurations with highest joint posterior probabilitiy among all sampled configurations (#regulators x #conditions).
     */
    virtual arma::mat& getMAP();
    
    virtual int getrand(int n);
	
};	 


#endif
