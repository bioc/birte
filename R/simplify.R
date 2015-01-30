simplify = function(affinities, cutoff=0.9){
  join.regulators = function(reg.type){
    reg = names(affinities[[reg.type]])
    genesets = sapply(affinities[[reg.type]], names)
    all.genes = unique(unlist(sapply(affinities[[reg.type]], names)))
    X = matrix(0, ncol=length(reg), nrow=length(all.genes))
    dimnames(X) = list(all.genes, reg)
    for(r in reg){
      X[genesets[[r]], r] = 1
    }
    D = dist(t(X), method="binary")
    hc = hclust(D, method="single")
    plot(hc)
    cl = cutree(hc, h=1-cutoff)
    maxcl = max(cl)
    cat("merging", length(reg), "regulators into", maxcl, "clusters\n")    
    a = lapply(1:maxcl, function(r){      
      g = unique(unlist(genesets[which(cl == r)]))
      v = rep(1, length(g))
      names(v) = g
      v
    })    
    names(a) = sapply(1:maxcl, function(r) paste(names(genesets)[which(cl == r)], sep="", collapse="_U_"))
    a
  }  
  
  reg.types = names(affinities)
  affinities = lapply(names(affinities)[sapply(affinities, length) > 0], join.regulators)
  names(affinities) = reg.types[sapply(affinities, length) > 0]
  affinities
}

proposeInteractions = function(affinities, cutoff.lower=0.1, cutoff.upper=0.8){  
  affinities = c(affinities$TF, affinities$miRNA)
  reg = names(affinities)
  genesets = sapply(affinities, names)
  all.genes = unique(unlist(sapply(affinities, names)))
  X = matrix(0, ncol=length(reg), nrow=length(all.genes))
  dimnames(X) = list(all.genes, reg)
  for(r in reg){
    X[genesets[[r]], r] = 1
  }
  Te = 1 - as.matrix(dist(t(X), method="binary"))
  Te[lower.tri(Te)] = 0
  take = which(Te > cutoff.lower & Te < cutoff.upper, arr.ind=TRUE)    
  cat("considering", NROW(take), "interaction terms\n")    
  a = lapply(1:NROW(take), function(i){
    int = intersect(genesets[[take[i,1]]], genesets[[take[i,2]]])
    v = rep(1, length(int))
    names(v) = int
    v
  })         
  names(a) = sapply(1:NROW(take), function(i) paste(names(genesets)[take[i,1]], "_", names(genesets)[take[i,2]], sep=""))  
  a
}
