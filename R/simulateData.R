# Simulate expression data
# 
# Author: frohlich



SigGenes <- function(rho, B = 100, m = 100) {
	SS = array(0, dim = c(B, B, m))
	for (h in 1:m) {
		if (h%%2 == 1) 
			r = rho
		else r = -rho
		SS[, , h] = toeplitz(r^(0:(B - 1)))
	}
	S <- diag(0, m * B)
	for (h in 1:m) for (i in 1:B) for (j in 1:B) 
				S[(h - 1) * B + i, (h - 1) * B + j] = SS[i, j, h]
	S
}


# simulate expression data for given states
sampleExpressionData = function(states=NULL, nrep=5, logFC=1){
	if(!is.null(states)){
		n = length(states)				
		logFC = logFC*states
	}
	else
		n = length(logFC)

	blocks = ceiling(n/100)
	SIG <- SigGenes(rho = 0.7, B = 100, m = blocks) / 10
	SIG = SIG[1:n, 1:n]				
	ex = t(mvrnorm(nrep, logFC, SIG))
	ex
}

sample.miRNAData = function(states, nrep=5, mean=0){	
	sampleExpressionData(states, nrep=nrep, logFC=mean)
}

sample.mRNAData = function(allgenes, states, types, genesets, nrep=5, fn.targets=0.1, fp.targets=0.2){		
	W = matrix(0, ncol=length(states), nrow=length(allgenes))
#   if(NCOL(interactions) > 0 && !is.null(interactions))
# 	  iterms = apply(interactions, 2, function(i) paste(i[1], "_", i[2], sep=""))
#   else
#     iterms = NULL
#	colnames(W) = c(names(genesets),iterms) 
  colnames(W) = names(genesets)
	rownames(W) = allgenes		
	#affinities = sapply(affinities, function(a) qnorm(1 - 10^(-abs(a)))) # convert into z-scores  
	cat("building W-matrix ...")	
	for(reg in names(genesets)){
		if(states[reg] == 1){
			regulon = setdiff(intersect(genesets[[reg]], allgenes), NA)						
			W[regulon, reg] = 1      
		}
	}	  
#   if(!is.null(interactions) && NCOL(interactions) > 0){
#   	for(i in 1:NCOL(interactions)){
#   		last = interactions[1,i]
#   		current = interactions[2,i]
#   		iterm = iterms[i]
#   		regulon1 = genesets[[last]]
#   		regulon2 = genesets[[current]]
#   		regulon = intersect(regulon1, regulon2)
#   		W[regulon, iterm] = 1
#   	}
# 	}
	# false positive target: no reaction, although expected
	if(fp.targets > 0){
		all.pred = which(W != 0)
		set0 = sample(all.pred, ceiling(length(all.pred)*fp.targets))
		W[set0] = 0
	}
	# false negative targets: reaction, although not expected
	if(fn.targets > 0){
		no.pred = setdiff(which(W == 0), set0)
		set1 = sample(no.pred, ceiling(length(no.pred)*fn.targets)) 
		W[set1] = 1 
	}
	mirs = which(types == "miRNA")
	W[, mirs] = -W[,mirs] # negative influence of miRNAs
	cat("done.\nBuilding expression matrix ...\n")
	# normalize: logFCs of -1 or 1 for differential genes expected
	W = t(apply(W, 1, function(x){
				if(sum(x) != 0)
					x/abs(sum(x[x!=0]))
				else
					x}))
#	states = c(states, rep(1, length(iterms)))	
	logFC = W %*% as.matrix(states)		
	control = sampleExpressionData(rep(0, length(logFC)), nrep=nrep, logFC=0)
	treated = sampleExpressionData(NULL, nrep=nrep, logFC=logFC)
	ex = cbind(control, treated)  
	print(head(ex))
	cat("done.\n")
	ex
}

# simulate data: general assumption is that in first condition (control) all miRNAs and TFs are inactive
simulateData = function(affinities, nrep=5, miRNAExpressions=TRUE, fn.targets=0.1, fp.targets=0.2, exp.nTF=5, exp.nmiR=5, exp.interact=5){		
	# sample miRNA and TF states for second condition	
	genesetsmiRNA = sapply(affinities$miRNA, names)	
	genesetsTFs = sapply(affinities$TF, names)	
  genesetsother = sapply(affinities$other, names)
	genesets = c(genesetsmiRNA, genesetsTFs)
	allgenes = unique(unlist(genesets))	
	interactions = miRNAstates = inter.states = NULL	
	
#	if(is.null(G)){
		repeat{
			miRNAstates = rbinom(length(genesetsmiRNA), 1, exp.nmiR/length(genesetsmiRNA))
			names(miRNAstates) =  names(genesetsmiRNA)
			if(sum(miRNAstates) > 0)
				break
		}
		# sample TF states for second condition
		repeat{
			TFstates = rbinom(length(genesetsTFs), 1, exp.nTF/length(genesetsTFs))
			names(TFstates) = names(genesetsTFs)
			if(sum(TFstates > 0))
				break
		}	
    # sample interaction effects
    if(!is.null(affinities$other)){      
      repeat{
        inter.states = rbinom(length(affinities$other), 1, exp.interact/length(affinities$other))
        names(inter.states) = names(affinities$other)
        if(sum(inter.states) > 0)
          break
      }
    }
    cat("==> simulated", sum(TFstates), "active TFs, ", sum(miRNAstates), "active miRNAs, ", sum(inter.states), "active interactions\n")
# 	}
# 	else{
# 		cat("random walk on graph ...")
# 		cl = clusters(G)
# 		G = induced.subgraph(G, V(G)$name[cl$membership == which.max(cl$csize)])
# 		current = sample(V(G)$name, 1)		
# 		TFstates = double(length(genesetsTFs))
# 		names(TFstates) = names(genesetsTFs)		
# 		miRNAstates = double(length(genesetsmiRNA))
# 		names(miRNAstates) = names(genesetsmiRNA)		
# 		or.logic = 1	
# 		last = current
# 		while(sum(TFstates) <= exp.nTF && sum(miRNAstates) <= exp.nmiR && NCOL(interactions) <= pmin(exp.nTF, exp.nmiR)){ # random walk of dem Graph
# 			if(or.logic == 1){
# 				if(current %in% names(miRNAstates))
# 					miRNAstates[current] = 1
# 				if(current %in% names(TFstates))
# 					TFstates[current] = 1
# 			}
# 			else{        
#  				m11 = match(last, interactions[1,])
#  				m12 = match(current, interactions[2,])
#  				m21 = match(last, interactions[2,])
#  				m22 = match(current, interactions[1,])
#  				if(all(is.na(c(m11, m12, m21, m22))) || (all(!is.na(c(m11, m12))) && m11 != m12) || (all(!is.na(c(m21, m22))) && m21 != m22))
#  					interactions = cbind(interactions, c(last, current))
#     	}				
# 			ad = unlist(neighborhood(G, 1, nodes=current))
# 			ad = setdiff(V(G)$name[ad], current)
# 			if(length(ad) == 0)
# 				current = sample(V(G)$name, 1)
# 			else{				
# 				last = current
# 				current = sample(ad, 1)
# 				regulon1 = genesets[[last]]
# 				regulon2 = genesets[[current]]
# 				conf.tab = rbind(cbind(length(intersect(regulon1, regulon2)), length(setdiff(regulon1, regulon2))), cbind(length(setdiff(regulon2, regulon1)), length(setdiff(allgenes, union(regulon1, regulon2)))))
# 				p = fisher.test(conf.tab, alternative="g")$p.value
# 				if(p < 0.05)
# 					or.logic = 0
# 				else
# 					or.logic = 1
# 			}
# 		}		
# 		cat("finished\n")
#     inter.states = double(NROW(interactions))
#     names(inter.sates) = paste(interactions[,1], interactions[,2], sep="_")
# 	}
	
	if(miRNAExpressions){
		# sample miRNA expression data	
		dat.miRNA = cbind(sample.miRNAData(miRNAstates, nrep=nrep, mean=0), sample.miRNAData(miRNAstates, nrep=nrep, mean=1))
		colnames(dat.miRNA) = c(rep("control", nrep), rep("treated", nrep))
		rownames(dat.miRNA) = names(genesetsmiRNA)    
	}
	else{
		dat.miRNA = NULL
	}	
	# sample TF expression data
	#s = sample(1:length(TFstates), round(state.err*length(TFstates)))	
  TFstates.obs = TFstates
	#TFstates.obs[s] = 1 - TFstates[s] # wrongly observed states due to low correlation with expression data
	dat.TF = cbind(sample.miRNAData(TFstates.obs, nrep=nrep, mean=0), sample.miRNAData(TFstates.obs, nrep=nrep, mean=0.5))	
	colnames(dat.TF) = c(rep("control", nrep), rep("treated", nrep))
	rownames(dat.TF) = names(genesetsTFs)
		
	# sample mRNA expression data2
	affinities = c(affinities$TF, affinities$miRNA, affinities$other)
	genesets = c(genesetsmiRNA, genesetsTFs, genesetsother)
  genesets = genesets[names(affinities)]
	states = c(miRNAstates, TFstates, inter.states)	
	states = states[match(names(genesets), names(states))]	  
  types = c(rep("miRNA", length(miRNAstates)), rep("TF", length(TFstates)), rep("interactions", length(inter.states)))
	dat.mRNA = sample.mRNAData(allgenes=allgenes, states=states, types=types, genesets=genesets, nrep=nrep, fn.targets=fn.targets, fp.targets=fp.targets)
	colnames(dat.mRNA) = c(rep("control", nrep), rep("treated", nrep))
	rownames(dat.mRNA) = allgenes	
  
	if(any(is.na(dat.mRNA)))
		stop("mRNA data contains NA!")
	list(dat.mRNA=dat.mRNA, dat.miRNA=dat.miRNA, dat.TF=dat.TF, miRNAstates=miRNAstates, TFstates=TFstates, inter.states=inter.states)
}

