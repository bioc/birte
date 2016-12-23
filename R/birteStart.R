
birteStart <- function(mRNAexpr, miRNAexpr=NULL, mRNA.data.type=c("array", "RNAseq"), miRNA.data.type=c("array", "RNAseq"), genesetsTF=NULL, genesetsmiRNA=NULL, genesetsothers=NULL, nrep.mRNA=c(5, 5), nrep.miRNA=NULL, nrep.TF=NULL, nrep.Q=NULL, alphamiR=1, betamiR=0.1, alpha_i=NULL, alpha_i0=NULL, niter=1e6, burnin=5*1e5, thin=50, model=c("all-plug-in", "no-plug-in"), only_switches=FALSE, noTF=FALSE, nomiRNA=FALSE, A_sigma=NULL, O_sigma, affinitiesTF=NULL, affinitiesmiRNA=NULL, affinitiesothers=NULL, potential_swaps=NULL, theta_TF=0.01, theta_miRNA=0.01, theta_other=0.01, K=NULL, init_S=NULL, init_T=NULL, init_other=NULL, TFexpr=NULL, TFexpr.data.type=c("array", "RNAseq"), alpha_i0TF=NULL, alpha_iTF=NULL, TF_sigma=NULL, alphaTF=NULL, betaTF=NULL, perc.overlap.cutoff=0.8, Qexpr=NULL, Qexpr.data.type=c("array", "RNAseq"), alpha_i0Q=NULL, alpha_iQ=NULL, Q_sigma=NULL, alphaQ=NULL, betaQ=NULL, alpha=NULL, beta=NULL, interactions=NULL) 
{
	model = match.arg(model, several.ok=FALSE)	
	model = switch(model, "all-plug-in"=1, "no-plug-in"=2)
	mRNA.data.type = match.arg(mRNA.data.type, several.ok=FALSE)
	miRNA.data.type = match.arg(miRNA.data.type, several.ok=FALSE)
	TFexpr.data.type = match.arg(TFexpr.data.type, several.ok=FALSE)
	Qexpr.data.type = match.arg(Qexpr.data.type, several.ok=FALSE)
	
	if(!is.null(miRNAexpr))
		miRNAexpr = as.matrix(miRNAexpr)
	if(!is.null(mRNAexpr))
		mRNAexpr = as.matrix(mRNAexpr)
	if(!is.null(Qexpr))
		Qexpr = as.matrix(Qexpr)
	
	C_cnt = length(nrep.mRNA)
	if(model == "all-plug-in"){
		if(!is.null(O_sigma) & length(O_sigma) != NROW(mRNAexpr))
			stop("length of O_sigma has to equal number of mRNAs!")
		if(!is.null(A_sigma) & length(A_sigma) != NROW(miRNAexpr))
			stop("length of A_sigma has to equal number of miRNAs!")
		if(!is.null(Q_sigma) & length(Q_sigma) != NROW(Qexpr))
			stop("length of Q_sigma has to equal number of other regulating factors!")
	}
	if(!is.null(alpha_i) & length(alpha_i) != NROW(miRNAexpr))
		stop("length of alpha_i has to equal number of miRNAs!")
	if(!is.null(alpha_i0) & length(alpha_i0) != NROW(miRNAexpr))
		stop("length of alpha_i0 has to equal number of miRNAs!")
	if(!is.null(alpha_iTF) & length(alpha_iTF) != NROW(TFexpr))
		stop("length of alpha_iTF has to equal number of TFs in TF data!")
	if(!is.null(alpha_i0TF) & length(alpha_i0TF) != NROW(TFexpr))
		stop("length of alpha_i0TF has to equal number of TFs in TF data!")
	if(!is.null(alpha_iQ) & length(alpha_iQ) != NROW(Qexpr))
		stop("length of alpha_iQ has to equal number of other regulating factors!")
	if(!is.null(alpha_i0Q) & length(alpha_i0Q) != NROW(Qexpr))
		stop("length of alpha_i0Q has to equal number of other regulating factors!")			
	if(!is.null(nrep.miRNA) && length(nrep.miRNA) != C_cnt)
		stop("Number of conditions provided for miRNA and mRNA data has to match!")		
	if(!is.null(nrep.TF) && length(nrep.TF) != C_cnt)
		stop("Number of conditions provided for TF and mRNA data has to match!")
	if(!is.null(nrep.Q) && length(nrep.Q) != C_cnt)
		stop("Number of conditions provided for other factors and mRNA data has to match!")		
	if(!is.null(init_T) && (NCOL(init_T) != length(genesetsTF) | NROW(init_T) != C_cnt))
		stop("init_T has to be of dimension #conditions x length of genesetsTF")
	if(!is.null(init_S) && (NCOL(init_S) != length(genesetsmiRNA) | NROW(init_S) != C_cnt))
		stop("init_S has to be of dimension #conditions x length of genesetsmiRNA")
	if(!is.null(init_other) && (NCOL(init_other) != length(genesetsothers) | NROW(init_other) != C_cnt))
		stop("init_other has to be of dimension #conditions x length of genesetsothers")
	if(!is.null(miRNAexpr) && length(genesetsmiRNA) != NROW(miRNAexpr)){
		stop("length of genesetsmiRNA and NROW(miRNAexpr) have to be equal!")
	}	
	if(!is.null(Qexpr) && length(genesetsothers) != NROW(Qexpr)){
		stop("length of genesetsothers and NROW(Qexpr) have to be equal!")
	}

	stopifnot(all(theta_miRNA > 0) & all(theta_TF > 0) & all(theta_other > 0) & alphamiR > 0 & betamiR > 0 & alphaTF > 0 & betaTF > 0 & alphaQ > 0 & betaQ > 0 & thin > 0 & burnin > 0 & niter > 0)
	
	miRNA = names(genesetsmiRNA)
	mRNA = rownames(mRNAexpr)
	TF = names(genesetsTF)	
	others = names(genesetsothers)
	A_cnt=as.integer(length(miRNA))
	T_cnt=as.integer(length(TF))
	others_cnt = as.integer(length(others))
	if(nomiRNA) {
		A_cnt = 0
		use_miRNA_expression = 0		
	}
	else{
		use_miRNA_expression = !is.null(miRNAexpr) & (NROW(miRNAexpr) > 0) & !is.null(nrep.miRNA)	
	}
	if(noTF) {
		T_cnt = 0
	}
	use_Q_expression = !is.null(Qexpr) & (NROW(Qexpr) > 0) & !is.null(nrep.Q)	
	
	## Edges from miRNAs to their targets
	mirTargets = genesetsmiRNA[miRNA]
	mirTargets = lapply(mirTargets, function(x) {which(mRNA %in% x)})
	## Edges from TFs to their targets
	TFtargets = genesetsTF[TF]
	TFtargets = lapply(TFtargets, function(x) {which(mRNA %in% x)})
	# Edges from others to their targets
	others.targets = genesetsothers[others]	
	others.targets = lapply(others.targets, function(x) {which(mRNA %in% x)})
	
	nTFexpr = 0
	if(! is.null(TFexpr)) {
		if(C_cnt > 1)
			nTFexpr = as.integer(dim(TFexpr)[1])
		else
			nTFexpr = as.integer(length(TFexpr))
	}			
	if(any(O_sigma == 0) | any(A_sigma == 0) | any(Q_sigma == 0))
		warning("Variance / dispersion parameter estimates contain 0s!")	
	if(!is.null(K)){
	  K_cnt = NCOL(K)
		K = as.numeric(K)		
	}
  else
    K_cnt = 0
			
  n0 = 1  
  if(length(theta_TF) < T_cnt && !is.null(theta_TF))
    theta_TF = rep(theta_TF[1], T_cnt)
	if(length(theta_miRNA) < A_cnt && !is.null(theta_miRNA))
	  theta_miRNA = rep(theta_miRNA[1], A_cnt)
	if(length(theta_other) < others_cnt && !is.null(theta_other))
	  theta_other = rep(theta_other[1], others_cnt)
	
	cat("\nBIRTE\n")
	cat("Data and network: #mRNAs = ", nrow(mRNAexpr), ", #miRNAs = ", A_cnt, ", #TFs = ", T_cnt, ", #other influence factors = ", others_cnt, ", #TFs with expression data = ", nTFexpr, ", miRNA expression data = ", use_miRNA_expression, ", use data for other influencing factors = ", use_Q_expression, "\n")
	if(C_cnt > 1){    
		#cat("Data type mRNA = ", mRNA.data.type, ", data type miRNA (if available) = ", miRNA.data.type, ", data type TFs (if available) = ", TFexpr.data.type, ", data type other regulators (if available) = ", Qexpr.data.type, "\n")    
		cat("#conditions = ", C_cnt, ", #replicates mRNA data = ", nrep.mRNA, "#replicates miRNA data = ", nrep.miRNA, "#replicates TF data = ", nrep.TF, ", #replicates other data = ", nrep.Q, "\n")
	}
	else{      
		cat("Using *relative* expression changes\n")	
	}
	cat("Prior parameters: theta_TF = ", head(theta_TF), ", theta_miRNA = ", head(theta_miRNA), ", theta_other = ", head(theta_other), "\n")	
  if(K_cnt > 0)
    cat("Using informative prior for regulator-regulator interactions\n")
	cat("Hyperparameters for variance prior (mRNA): alpha = ", alpha, ", beta = ", beta, "\n")
	if(model == 2 && (!is.null(miRNAexpr) || !is.null(TFexpr)))
		cat("Hyperparameters for variance priors (miRNA, TFs): alpha (miRNA) = ", alphamiR, ", beta (miRNA) = ", betamiR, ", alpha (TF) = ", alphaTF, ", beta (TF) = ", betaTF, "\n")
  if(model == 2 && !is.null(Qexpr))
    cat("Hyperparameters for variance prior (other data): alpha (other) = ", alphaQ, ", beta (other) = ", betaQ, "\n")
	cat("MCMC parameters: burnin = ", burnin, ", niter = ", niter, ", thin = ", thin, ", MODEL = ", model, "\n\n")	  
	replicates = c(nrep.miRNA, nrep.mRNA, nrep.TF, nrep.Q)		    
	result = .Call("getStates", num_conditions=as.integer(C_cnt), nmRNA=as.integer(length(mRNA)), nmiRNA=as.integer(A_cnt), nTF=as.integer(T_cnt), nothers=as.integer(others_cnt), replicates=as.integer(replicates), 
					  mRNA_expression=as.numeric(mRNAexpr), miRNA_expression=as.numeric(miRNAexpr), TFexpr=as.numeric(TFexpr), Qexpr=as.numeric(Qexpr), nTFexpr=as.integer(nTFexpr),
					  O_sigma=as.numeric(O_sigma), A_sigma=as.numeric(A_sigma), TF_sigma=as.numeric(TF_sigma), Q_sigma=as.numeric(Q_sigma),
					  mRNADataType=as.integer((mRNA.data.type=="RNAseq")*1), miRNADataType=as.integer((miRNA.data.type=="RNAseq")*1), TFexprDataType=as.integer((TFexpr.data.type=="RNAseq")*1), QexprDataType=as.integer((Qexpr.data.type=="RNAseq")*1),
				    use_miRNA_expression=as.integer(use_miRNA_expression), use_Q_expression=as.integer(use_Q_expression),
					  
					  genesetsmiRNA=mirTargets, genesetsTF=TFtargets, genesetsothers=others.targets,
					  alpha=as.numeric(alpha), beta=as.numeric(beta),
				    n0 = as.numeric(n0), alphamiR=as.numeric(alphamiR), betamiR=as.numeric(betamiR), alphaTF=as.numeric(alphaTF), betaTF=as.numeric(betaTF), alphaQ=as.numeric(alphaQ), betaQ=as.numeric(betaQ),
				    alpha_i0 = as.numeric(alpha_i0), alpha_i = as.numeric(alpha_i),  alpha_i0TF=as.numeric(alpha_i0TF), alpha_iTF=as.numeric(alpha_iTF),  alpha_i0Q=as.numeric(alpha_i0Q), alpha_iQ=as.numeric(alpha_iQ), 				       					 
				      
					  model=as.integer(model), niter=as.integer(niter), burnin=as.integer(burnin), thin=as.integer(thin), 
					  
					  only_switches=as.integer(only_switches), T_potential_swaps = sapply(potential_swaps$TF, as.integer), S_potential_swaps = sapply(potential_swaps$miRNA, as.integer), Q_potential_swaps=sapply(potential_swaps$other, as.integer), 					 
					  theta_TF=as.numeric(theta_TF), theta_miRNA=as.numeric(theta_miRNA), theta_other=as.numeric(theta_other), K=K, K_cnt=as.integer(K_cnt), interactions=as.integer(interactions), 
					  init_S=as.integer(init_S), init_T=as.integer(init_T), init_other=as.integer(init_other),
					  
					  affinitiesTF=affinitiesTF, PTC_miRNA=affinitiesmiRNA, affinitiesothers=affinitiesothers, PACKAGE="birte")		

  if(!is(result$post, "matrix")){
	  result$post = matrix(result$post, ncol=C_cnt)    
  }
	if(!is(result$map, "matrix")){
	  result$map = matrix(result$map, ncol=C_cnt)    
	}
  dimnames(result$post) = list(c(names(affinitiesTF), names(affinitiesmiRNA), names(affinitiesothers)), paste("condition", 1:C_cnt))
  dimnames(result$map) = dimnames(result$post)
	result$contains.interactions = !is.null(affinitiesothers) & !is.null(K)  
  affinities = c(affinitiesTF, affinitiesmiRNA, affinitiesothers)
	result$design = matrix(0, ncol=length(affinities) + 1, nrow=NROW(mRNAexpr))
  dimnames(result$design) = list(rownames(mRNAexpr), c("intercept", names(affinities)))
  result$design[,1] = 1
  for(i in setdiff(colnames(result$design), "intercept")){
    result$design[intersect(rownames(mRNAexpr), names(affinities[[i]])), i] = affinities[[i]][intersect(rownames(mRNAexpr), names(affinities[[i]]))]
  } 	
  for(c in 1:C_cnt){
    for(r in 1:NCOL(result$coef)){
      if(NROW(result$coef[c,r][[1]]) > 0)
        rownames(result$coef[c,r][[1]]) = colnames(result$design)    
    }
  }  	
	result$nburnin = burnin  
  result$C_cnt = C_cnt;  
  result$param = list(alpha=alpha, beta=beta)  
	return(result)
}
	
