# Main user interface functions
# Author: Holger Fröhlich
###############################################################################

birte.limma.aux = function(limmamiRNA, dat.miRNA, lfc, fdr, explain.LFC=TRUE){
	if(!is.null(limmamiRNA) && !is.null(dat.miRNA)){
		cat("Automatic assignement of #replicates with design matrix.\n")
		stopifnot(!is.null(limmamiRNA$design))	
		if(is.null(limmamiRNA$pvalue.tab$ID))
		  limmamiRNA$pvalue.tab$ID = rownames(limmamiRNA$pvalue.tab)
    if(!explain.LFC){
      stopifnot(NCOL(limmamiRNA$design) == 2)
		  groups = lapply(1:NCOL(limmamiRNA$design), function(i) which(limmamiRNA$design[,i] == 1))		
		  dat.miRNA[,as.vector(unlist(groups))]
		  nrep.miRNA = sapply(groups, length)		
    }    		 
    else{
      dat.miRNA = as.matrix(limmamiRNA$pvalue.tab$logFC)
      rownames(dat.miRNA) = limmamiRNA$pvalue.tab$ID
      nrep.miRNA = c(1,1)
    }
		stopifnot(!is.null(dat.miRNA))
		if(is.null(lfc))
			lfc = 0
		if(is.null(fdr))
			fdr = 0.05   
		diff.miRNA = as.character(limmamiRNA$pvalue.tab$ID[abs(limmamiRNA$pvalue.tab$logFC) > lfc & limmamiRNA$pvalue.tab$adj.P.Val < fdr])
    
		A_Sigma = limmamiRNA$lm.fit$sigma[rownames(dat.miRNA)]    
	}
	else{
		cat("No data available (limma object is NULL).\n")
		dat.miRNA = NULL
	}
	list(nrep=nrep.miRNA, sigma=A_Sigma, dat=dat.miRNA, diff.exp=diff.miRNA)
}

# convenience function using limma
birteLimma = function(dat.mRNA, limmamRNA,
		data.regulators=NULL, limma.regulators=NULL, fdr.regulators=NULL, lfc.regulators=NULL, 
		init.regulators=NULL, theta.regulators=NULL, reg.interactions=FALSE, affinities, use.affinities=FALSE,
		niter=100000, nburnin=100000, thin=50, potential_swaps=NULL, only_switches=FALSE, only.diff.TFs=TRUE, explain.LFC=TRUE, model=c("no-plug-in", "all-plug-in")){

  if(is.null(limmamRNA))
    stop("Please provide argument limmamRNA!")
  cat("Automatic assignement of #replicates with design matrix (mRNA data).\n")
  stopifnot(!is.null(limmamRNA$design))
  if(is.null(limmamRNA$pvalue.tab$ID))
    limmamRNA$pvalue.tab$ID = rownames(limmamRNA$pvalue.tab)
  if(!explain.LFC){
    stopifnot(NCOL(limmamRNA$design) == 2)
    groups = lapply(1:NCOL(limmamRNA$design), function(i) which(limmamRNA$design[,i] == 1))  
    dat.mRNA = dat.mRNA[,as.vector(unlist(groups))]
    nrep.mRNA = sapply(groups, length)	
  }
  else{
    dat.mRNA = as.matrix(limmamRNA$pvalue.tab$logFC)
    rownames(dat.mRNA) = limmamRNA$pvalue.tab$ID
    nrep.mRNA = c(1,1)
  }
  stopifnot(!is.null(dat.mRNA))		
  O_Sigma = limmamRNA$lm.fit$sigma[rownames(dat.mRNA)]	
  df.mRNA = limmamRNA$lm.fit$df.residual
  names(O_Sigma) = names(df.mRNA) = rownames(dat.mRNA)	
  
  dat = diff.exp = sigma = nrep = list()
  for(regulator.type in names(limma.regulators)){
    cat(regulator.type, "==> ")
    if(!is.null(limma.regulators[[regulator.type]])){
      res = birte.limma.aux(limma.regulators[[regulator.type]], data.regulators[[regulator.type]], lfc.regulators[[regulator.type]], fdr.regulators[[regulator.type]], explain.LFC=explain.LFC)
    }
    else
      res = NULL
    dat[[regulator.type]] = res$dat
    sigma[[regulator.type]] = res$sigma    
    nrep[[regulator.type]] = res$nrep
    diff.exp[[regulator.type]] = res$diff.exp  	
  }	
  
	res = birteRun(dat.mRNA=dat.mRNA, mRNA.Sigma=O_Sigma, nrep.mRNA=nrep.mRNA, df.mRNA=df.mRNA,
					data.regulators=dat, sigma.regulators=sigma,nrep.regulators=nrep, diff.regulators = diff.exp,
					init.regulators=init.regulators, theta.regulators=theta.regulators,	reg.interactions=reg.interactions, 				
					affinities=affinities, use.affinities=use.affinities, 
					niter=niter, nburnin=nburnin, thin=thin, potential_swaps=potential_swaps, only_switches=only_switches, only.diff.TFs=only.diff.TFs, explain.LFC=explain.LFC, model=model)
	res
}

# main function
birteRun = function(dat.mRNA, mRNA.Sigma=NULL, nrep.mRNA=c(5, 5), df.mRNA=sum(nrep.mRNA)-2,
		data.regulators=NULL, sigma.regulators=NULL, nrep.regulators=NULL, diff.regulators=NULL, 
		init.regulators=NULL, theta.regulators=NULL, reg.interactions=FALSE, affinities, use.affinities=FALSE,
    niter=100000, nburnin=100000, thin=50, potential_swaps=NULL, only_switches=FALSE, only.diff.TFs=TRUE, explain.LFC=TRUE, model=c("no-plug-in", "all-plug-in")) {
  
  
  model = match.arg(model, several.ok=FALSE)  
	if(is.null(affinities))
		stop("Please provide regulator binding affinities! The TF-target and miRNA-target networks are required to run biRte!")	
	if(length(affinities) > 3)
		stop("biRte right now only supports up to 3 regulator types and corresponding datasets!")
	ndat = length(data.regulators)		
  nreg = length(affinities)
	if(length(nrep.regulators)!=ndat || length(diff.regulators) != ndat ||  (!is.null(sigma.regulators) && length(sigma.regulators) != ndat) || (!is.null(init.regulators) && length(init.regulators) != nreg))
		stop("Please provide for each regulator dataset at least number of replicates per condition and differentially expressed regulators!\nIf initial states for regulators are provided these initializations need to be done for all regulator types!")
	#if(ndat > 0 && length(affinities) != ndat)
	#	stop("Number of components of affinities have to match number of components of data.regulators!")
	if(!(any(c("TF", "miRNA", "other") %in% names(affinities))))
		stop("Regulator types have to be one of 'TF', 'miRNA', 'other'")
	stopifnot(all(names(data.regulators) %in% names(affinities)))
  if(length(nrep.mRNA) != 2)
    stop("biRte assumes exactly two conditions!")
	#
	cat("Formatting regulator-target network -> checking overlap between network and measurements.\n")	
	C_cnt = length(nrep.mRNA)
	alpha0 = alpha = alpha.para = beta.para = genesets = theta.regulators2 = init.regulators2 = data.regulators2 = regulators.data.type = potential_swaps2 = sigma.regulators2 = logFC.reg = list()
	alltargets = c()
	all.regulators = c()	
  affinities.orig = affinities
  if(explain.LFC)
    cat("--> biRte tries to explain mRNA log fold changes\n")
  else
    cat("--> biRte tries to explain condition specific mRNA expression\n")
  if(is.null(mRNA.Sigma) && all(nrep.mRNA > 1) && NCOL(dat.mRNA) > 1){
    warning("No SDs for mRNA data provided --> trying to automatically build limma model and deduce standard deviations ...")
    colnames(dat.mRNA) = c(rep("condition1", nrep.mRNA[1]), rep("condition2", nrep.mRNA[2]))    
    limma.mRNA = limmaAnalysis(dat.mRNA, contrast="condition1 - condition2")
    mRNA.Sigma = limma.mRNA$lm.fit$sigma        
    names(mRNA.Sigma) = rownames(dat.mRNA)
    if(any(limma.mRNA$lm.fit$df.residual != df.mRNA)){
      df.mRNA = limma.mRNA$lm.fit$df.residual
      names(df.mRNA) = names(mRNA.Sigma)
      warning("Degrees of freedom provided for mRNA data don't agree with those deduced from limma model! Please check your inputs!")
    }      
  }
  if(length(df.mRNA) == 1){
    df.mRNA = rep(df.mRNA, NROW(dat.mRNA))    
    names(df.mRNA) = rownames(dat.mRNA)    
  }
	for(regulator.type in names(affinities)){		
    cat(regulator.type, "\n")    
		myaffinities = affinities[[regulator.type]]		
		res = birte.aux(model=model, regulator.type=regulator.type, C_cnt=C_cnt, nrep.mRNA=nrep.mRNA, dat.mRNA=dat.mRNA, dat.miRNA=data.regulators[[regulator.type]], 
				miRNA.data.type="array", miRNA.Sigma=sigma.regulators[[regulator.type]], nrep.miRNA=nrep.regulators[[regulator.type]], diff.miRNA=diff.regulators[[regulator.type]], 
				init_miR=init.regulators[[regulator.type]], theta_miRNA=theta.regulators[[regulator.type]], affinitiesmiRNA=myaffinities, use.affinities=use.affinities, only.switches=only_switches, potential_swaps=potential_swaps, only.diff.TFs=only.diff.TFs, explain.LFC=explain.LFC)
		data.regulators2[[regulator.type]] = res$dat		
		alpha0[[regulator.type]] = res$alpha0
		alpha[[regulator.type]] = res$alpha
		alpha.para[[regulator.type]] = res$alpha.para
		beta.para[[regulator.type]] = res$beta.para
		affinities[[regulator.type]] = res$affinities        
		genesets[[regulator.type]] = sapply(res$affinities, names)
		regulators.data.type[[regulator.type]] = "array"
		theta.regulators2[[regulator.type]] = res$theta
		init.regulators2[[regulator.type]] = res$init
		potential_swaps2[[regulator.type]] = res$potential_swaps
		sigma.regulators2[[regulator.type]] = res$sigma
    logFC.reg[[regulator.type]] = res$logFC    
		if(explain.LFC){
			if(!is.null(res$dat)){
				nrep = nrep.regulators[[regulator.type]][1]
				n = NCOL(data.regulators[[regulator.type]])
				data.regulators2[[regulator.type]] = logFC.reg[[regulator.type]]
			}
			nrep.regulators[[regulator.type]] = 1
			init.regulators2[[regulator.type]] = init.regulators2[[regulator.type]][1,,drop=FALSE] 
      alpha0[[regulator.type]] = res$alpha0*0      
		}
		#		
		alltargets = union(alltargets, unlist(genesets[regulator.type]))    
		ctrl = setdiff(unlist(genesets[[regulator.type]]), rownames(dat.mRNA))    
		if(length(ctrl) > 0)
			stop(paste("Problem: Incompatibilities between network and annnotation of", regulator.type, "data (not the same target genes)"))
		all.regulators = union(all.regulators, colnames(init.regulators2[[regulator.type]]))
	}
  #  
  reg.type.nprov = setdiff(c("TF", "miRNA", "other"), names(affinities))
  for(r in reg.type.nprov){
    if(explain.LFC)
      nrep.regulators[[r]] = 0
    else	
      nrep.regulators[[r]] = rep(0, C_cnt)
  }
  nrep.regulators = as.list(nrep.regulators)	
  #
  K = NULL	
	interactions=NULL  
  if(!is.null(affinities$other) && reg.interactions){
    if(!is.null(theta.regulators$other) && is(theta.regulators$other, "matrix")){
      K = theta.regulators$other
      theta.regulators$other = NULL
    }
    com = strsplit(names(affinities$other), "_")      	  
    regulators.interact = unique(unlist(com))
    if(!is.null(K))
      K = K[intersect(regulators.interact, rownames(K)), intersect(regulators.interact, colnames(K))]
    else{
      K = matrix(1/choose(length(regulators.interact),2), nrow=length(regulators.interact), ncol=length(regulators.interact))
      dimnames(K) = list(regulators.interact, regulators.interact)      
    }			
    interactions = t(sapply(com, function(x) c(which(rownames(K) == x[1]), which(rownames(K) == x[2])))) - 1        	    
  }    
  #  
	dat.mRNA = dat.mRNA[alltargets,,drop=FALSE]		    
	mRNA.Sigma = mRNA.Sigma[alltargets]
  df.mRNA = df.mRNA[alltargets]
	if(any(is.na(mRNA.Sigma)))
		stop("Variance of mRNA expressions must not equal NA!")
	stopifnot(rownames(dat.mRNA) == rownames(mRNA.Sigma))	
  
  dat.mRNA.orig = dat.mRNA  
  est = suppressWarnings(fitdistr(1/mRNA.Sigma^2, "gamma"))
	alpha.mRNA = est$estimate[1]
	beta.mRNA = est$estimate[2]    
	var.post = (2*beta.mRNA + df.mRNA*mRNA.Sigma^2) / (2*alpha.mRNA + df.mRNA) # limma estimate of posterior mean
	if(explain.LFC){
    if(all(nrep.mRNA == 1)) # log FCs provided
      logFC = dat.mRNA
    else{
	    mu = cbind(apply(dat.mRNA[,1:nrep.mRNA[1]], 1, mean), apply(dat.mRNA[,(nrep.mRNA[1] + 1):NCOL(dat.mRNA)], 1, mean))
	    rownames(mu) = rownames(dat.mRNA)  	
		  logFC = mu[,2] - mu[,1]	    
    }
		dat.mRNA = logFC / var.post
		nrep.mRNA = 1					    
	}	
  else{
    dat.mRNA = t(sapply(1:NROW(dat.mRNA), function(i) dat.mRNA[i,] / var.post[i]))
  }	
  # post normalization:
  beta.mRNA = alpha.mRNA # expectation=1 and variance = 1/alpha
	res = birteStart(mRNAexpr=dat.mRNA, miRNAexpr=data.regulators2$miRNA, Qexpr=data.regulators2$other, nrep.mRNA=nrep.mRNA, nrep.miRNA=nrep.regulators$miRNA, nrep.TF=nrep.regulators$TF, nrep.Q=nrep.regulators$other,
			mRNA.data.type="array", miRNA.data.type=regulators.data.type$miRNA, Qexpr.data.type=regulators.data.type$other,
			genesetsTF=genesets$TF, genesetsmiRNA=genesets$miRNA, genesetsothers=genesets$other,
			alpha_i=alpha$miRNA, alpha_i0=alpha0$miRNA, alpha_i0Q=alpha0$other, alpha_iQ=alpha$other,
			alphamiR=alpha.para$miRNA, betamiR=beta.para$miRNA, alphaQ=alpha.para$other, betaQ=beta.para$other,
			niter=niter, burnin=nburnin, thin=thin, model=model, 
			only_switches=only_switches, nomiRNA=is.null(genesets$miRNA), noTF=is.null(genesets$TF), affinitiesTF=affinities$TF, affinitiesmiRNA=affinities$miRNA, affinitiesothers=affinities$other,
			potential_swaps=potential_swaps2, 
			theta_TF=theta.regulators2$TF, theta_miRNA=theta.regulators2$miRNA, theta_other=theta.regulators2$other, K=K, interactions=interactions,
			A_sigma=sigma.regulators2$miRNA, O_sigma=mRNA.Sigma, Q_sigma=sigma.regulators2$other, init_S=init.regulators2$miRNA, init_T=init.regulators2$TF, init_other=init.regulators2$other, 
			TFexpr=data.regulators2$TF, TFexpr.data.type=regulators.data.type$TF, 
			alpha_i0TF=alpha0$TF, alpha_iTF=alpha$TF, TF_sigma=sigma.regulators2$TF, alphaTF=alpha.para$TF, betaTF=beta.para$TF, alpha=alpha.mRNA, beta=beta.mRNA)	
	#
  res$affinities = affinities.orig  	
  res$explain.LFC = explain.LFC  
  res$df.mRNA = df.mRNA
  res
}

#  auxilliary function
birte.aux = function(model="all-plug-in", regulator.type, C_cnt, dat.mRNA, nrep.mRNA, dat.miRNA, miRNA.data.type="array", miRNA.Sigma, nrep.miRNA, diff.miRNA, init_miR, theta_miRNA, affinitiesmiRNA, use.affinities, only.switches, potential_swaps, only.diff.TFs, explain.LFC){			
	alpha_i0 = alpha.miRNA = alpha.miR = beta.miR = miRNA.data.type = NULL	
	if(length(affinitiesmiRNA) > 0 ) {
	  affinitiesmiRNA = sapply(affinitiesmiRNA, function(s) s[intersect(names(s), rownames(dat.mRNA))])		    
		if(length(affinitiesmiRNA) == 0)
			warning(paste("No overlap between targets of type", regulator.type, "and mRNA measurements!"))
		
		if(any(sapply(affinitiesmiRNA, length) == 0)) {
			warning("Not all regulator genesets have non-zero length --> removing empty genesets (check result$affinities)!")
			affinitiesmiRNA = affinitiesmiRNA[sapply(affinitiesmiRNA, length) > 0]			
		}	
    if(!use.affinities){
      if(regulator.type == "miRNA")
  		  affinitiesmiRNA = lapply(affinitiesmiRNA, function(a) -a/a) 
      else
        affinitiesmiRNA = lapply(affinitiesmiRNA, function(a) a/a) 
    }
    
		if(is.null(init_miR)){
			init_miR = matrix(0, nrow=C_cnt, ncol=length(affinitiesmiRNA))
			colnames(init_miR) = names(affinitiesmiRNA)
      init_miR[intersect(diff.miRNA, names(init_miR))] = 1
		}
    if(!is(init_miR, "matrix")){
      init_miR2 = matrix(0,  nrow=C_cnt, ncol=length(affinitiesmiRNA))
      colnames(init_miR2) = names(affinitiesmiRNA)
      init_miR2[1,] = init_miR[names(affinitiesmiRNA)]
      init_miR = init_miR2
    }
		init_miR = init_miR[,names(affinitiesmiRNA), drop=FALSE]		
    logFC = NULL    
		if(!is.null(dat.miRNA)){			
			miRNA.data.type = match.arg(miRNA.data.type, several.ok=FALSE)
			miRNAInAnnotation = intersect(names(affinitiesmiRNA), rownames(dat.miRNA))
			dat.miRNA = dat.miRNA[miRNAInAnnotation,,drop=FALSE]							
			miRNA.Sigma = miRNA.Sigma[miRNAInAnnotation]
			if(regulator.type != "TF"){				
				affinitiesmiRNA = affinitiesmiRNA[miRNAInAnnotation]
				init_miR = init_miR[,miRNAInAnnotation, drop=FALSE]
			}          			
			if(is.null(miRNA.Sigma)){        
			  warning("No SDs for regulator type ", regulator.type, " provided --> trying to automatically build limma model and deduce standard deviations ...")
        if(all(nrep.miRNA == 1))
          stop("Cannot build limma model ==> #replicates = 1!")
			  colnames(dat.miRNA) = c(rep("condition1", nrep.miRNA[1]), rep("condition2", nrep.miRNA[2]))    
			  limma.miRNA = limmaAnalysis(dat.miRNA, contrast="condition2 - condition1")
			  miRNA.Sigma = limma.miRNA$lm.fit$sigma
			}
			
			alpha_i0 = apply(dat.miRNA[,1:nrep.miRNA[1], drop=FALSE], 1, mean, na.rm=TRUE) # average expression in first condition
			names(alpha_i0) = rownames(dat.miRNA)
			if(model == "all-plug-in")
				stopifnot(rownames(dat.miRNA) == rownames(miRNA.Sigma))
			stopifnot(length(nrep.miRNA) == length(nrep.mRNA))			
			cat(length(diff.miRNA), " differential regulators of type" , regulator.type, "found\n")	     	
			if(length(nrep.miRNA) != 2){
				warning(paste(regulator.type, "data does not contain 2 conditions! It will be ignored for model inference!\n"))
				dat.miRNA = NULL				
			}
			else{				
        if(all(nrep.miRNA == 1) && explain.LFC) # logFC data provided
          logFC = dat.miRNA
        else{
			    conditions = as.vector(unlist(sapply(1:length(nrep.miRNA), function(co) rep(paste("condition", co, sep=""), nrep.miRNA[co]))))  	
				  logFC = as.matrix(apply(dat.miRNA[, conditions == "condition2", drop=FALSE], 1, mean) - apply(dat.miRNA[, conditions == "condition1", drop=FALSE], 1, mean))
          rownames(logFC) = rownames(dat.miRNA)
        }        
				est = suppressWarnings(fitdistr(1/miRNA.Sigma^2, "gamma")$estimate)
				alpha.miR = est[1]
				beta.miR = est[2]
				which.up = rownames(logFC)[logFC > 0]
				which.down = rownames(logFC)[logFC < 0]				
				if(length(diff.miRNA) > 0){
					stopifnot(class(diff.miRNA) == "character")		
					diff.miRNA = intersect(diff.miRNA, rownames(dat.miRNA))
					diff.miRNA.up = intersect(rownames(logFC)[logFC > 0], diff.miRNA)
					diff.miRNA.down = intersect(rownames(logFC)[logFC < 0], diff.miRNA)
				}
				else{
					diff.miRNA.up = diff.miRNA.down = c()
				}
				alpha.miRNA = rep(1, nrow(dat.miRNA))
				names(alpha.miRNA) = rownames(dat.miRNA)
				alpha.miRNA[which.down] = -1
				if(length(diff.miRNA.up) > 0){																
					logFC.up = drop(logFC[diff.miRNA.up,])
					alpha.miRNA[which.up] = median(logFC.up, na.rm=TRUE)
					alpha.miRNA[diff.miRNA.up] = logFC.up
				}
				if(length(diff.miRNA.down) > 0){																
					logFC.down = drop(logFC[diff.miRNA.down,])
					alpha.miRNA[which.down] = median(logFC.down, na.rm=TRUE)
					alpha.miRNA[diff.miRNA.down] = logFC.down
				}					
				if(length(diff.miRNA.up) + length(diff.miRNA.down) > 0 && regulator.type == "TF" && only.diff.TFs){ # für TFs: nur differentielle werden berücksichtigt!
					dat.miRNA = dat.miRNA[union(diff.miRNA.up, diff.miRNA.down),,drop=FALSE]
          logFC = logFC[rownames(dat.miRNA),,drop=FALSE]					
					miRNA.Sigma = miRNA.Sigma[rownames(dat.miRNA)]
					alpha.miRNA = alpha.miRNA[rownames(dat.miRNA)]
					alpha_i0 = alpha_i0[rownames(dat.miRNA)]
          affinitiesmiRNA = affinitiesmiRNA[c(rownames(dat.miRNA), setdiff(names(affinitiesmiRNA), rownames(dat.miRNA)))] # rotiere Regulatoren mit Daten an den Anfang!!!          
					init_miR = init_miR[,names(affinitiesmiRNA), drop=FALSE]                                    
				}				
				stopifnot(all(miRNA.Sigma != 0))
			}							
		}	
		else
			nrep.miRNA = rep(0, C_cnt)
		if(is.null(theta_miRNA)){
			theta_miRNA = rep(1/length(affinitiesmiRNA), length(affinitiesmiRNA))
      names(theta_miRNA) = names(affinitiesmiRNA)
      theta_miRNA[intersect(diff.miRNA, names(theta_miRNA))] = 0.8
		}
	}
	else{
		nrep.miRNA = rep(0, C_cnt)
		theta_miRNA = 1e-10
	}	
	if(length(theta_miRNA) > 1)
	  theta_miRNA = theta_miRNA[names(affinitiesmiRNA)]  
	if(!only.switches && (is.null(potential_swaps) || length(potential_swaps[[regulator.type]]) != length(affinitiesmiRNA))){
	  cat("Calculating potential swaps for regulator type", regulator.type, "\n")
	  mypotential_swaps = getPotentialSwaps(sapply(affinitiesmiRNA, names))
	}
	else
	  mypotential_swaps = potential_swaps[[regulator.type]][names(affinitiesmiRNA)]
  
	list(dat=dat.miRNA, sigma=miRNA.Sigma, logFC=logFC, data.type=miRNA.data.type, nrep=nrep.miRNA, init=init_miR, theta=theta_miRNA, alpha0=alpha_i0, alpha=alpha.miRNA, alpha.para=alpha.miR, beta.para=beta.miR, affinities=affinitiesmiRNA, potential_swaps=mypotential_swaps)
}
