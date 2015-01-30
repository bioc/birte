estimateNetwork = function(model, thresh=0.1, select=c("marginal", "MAP"), method="pairwise", de.genes, bootstrap=0, typeII=0.1){    
  if(model$C_cnt > 2)
    stop("Method does only work for two conditions!")    
  stopifnot(!is.null(de.genes))
  select = match.arg(select, several.ok=FALSE)  
  affinities = c(sapply(model$affinities$TF, names), sapply(model$affinities$miRNA, names), sapply(model$affinities$other, names))
  regulators = list()
  if(select == "marginal"){
    for(c in 1:model$C_cnt){
      regulators[[c]] = rownames(model$post)[model$post[,c] > thresh]          
    }
  }
  else{
    for(c in 1:model$C_cnt){
      regulators[[c]] = rownames(model$map)[model$map[,c] > 0]          
    }
  }
  if(!model$explain.LFC)
    regulators[[1]] = setdiff(regulators[[1]], regulators[[2]])      
  myregulators = setdiff(unlist(strsplit(regulators[[1]], "_")), "U")
  myregulators = myregulators[myregulators %in% names(affinities)]  
  alltargets = intersect(unlist(affinities[myregulators]), de.genes)    
  fp.targets = sapply(affinities[myregulators], function(x) floor(0.05*length(intersect(de.genes, x))))
  X = matrix(0, nrow=length(alltargets), ncol=length(myregulators))
  colnames(X) = myregulators
  rownames(X) = alltargets   
  for(reg in myregulators){
    X[intersect(affinities[[reg]], rownames(X)), reg] = 1    
  }        
  cs = colSums(X)  
  X = X[, cs > fp.targets, drop=FALSE]
  nem.model = NULL
  if(NCOL(X) > 1){
    control = set.default.parameters(colnames(X), para=c(0.05, typeII))      
    if(bootstrap == 0)
        nem.model = nem(X, inference=method, control=control)      
    else
      nem.model = nem.bootstrap(X, inference=method, control=control, nboot=bootstrap)
  }
  else
    cat("Insufficient number of differentially expressed target genes ==> network inference impossible.")
  
  nem.model
}
