birtePredict = function(model, test.genes, method=c("Bayes", "MAP"), knock.out=NULL){  
  method = match.arg(method)
  x.test = matrix(0, ncol=ncol(model$design), nrow=length(test.genes))
  dimnames(x.test) = list(test.genes, colnames(model$design))
  affinities = c(sapply(model$affinities$TF, names), sapply(model$affinities$miRNA, names), sapply(model$affinities$other, names))
  types = c(rep("TF", length(model$affinities$TF)), rep("miRNA", length(model$affinities$miRNA)), rep("other", length(model$affinities$other)))
  names(types) = names(affinities)
  x.test[,1] = 1
  for(i in setdiff(colnames(x.test), union("intercept", knock.out))){
    if(types[i] == "miRNA")
      x.test[intersect(test.genes, affinities[[i]]), i] = -1
    else
      x.test[intersect(test.genes, affinities[[i]]), i] = 1
  } 
  all.coef = c("intercept", rownames(model$post))
  pred = list()
  for(c in 1:model$C_cnt){    
    pred[[c]] = list()    
    if(method == "Bayes"){
      for(r in 1:NCOL(model$coef)){       
        predtmp = x.test[,all.coef] %*% (model$coef[c,r][[1]][all.coef,] * c(1, model$post[,c]))        
        pred[[c]][[r]] = data.frame(gene=test.genes, mean=apply(predtmp, 1, mean), sd=apply(predtmp, 1, sd)) 
      }
    }
    else if(method == "MAP"){
      if(!is.null(model$fit.ridge)){            
        if(length(setdiff(names(coef(model$fit.ridge)), "(Intercept)")) > 0)
          predtmp = predict(model$fit.ridge, data.frame(x.test[, intersect(colnames(x.test), names(coef(model$fit.ridge)))]))
        else
          predtmp = predict(model$fit.ridge)
        pred[[c]] = data.frame(gene=test.genes, mean=predtmp)      
      }
      else
        stop("No predictions possible: either use Bayesian prediction or fit ridge regression model first!")
    }          
  }
  pred
}

birteFitRidge = function(model, mRNA.train, ref.cond=1){
  affinities = c(sapply(model$affinities$TF, names), sapply(model$affinities$miRNA, names), sapply(model$affinities$other, names))
  map = model$map[, ref.cond]
  map = map[map != 0]  
  if(class(mRNA.train) != "numeric")
    stop("Model can only fit one dimensional response variables!")  
  if(length(map) > 0){
    x.train = matrix(0, ncol=length(map), nrow=length(mRNA.train))
    dimnames(x.train) = list(names(mRNA.train), names(map))
    affinities = c(sapply(model$affinities$TF, names), sapply(model$affinities$miRNA, names), sapply(model$affinities$other, names))  
    for(i in colnames(x.train)){
      x.train[intersect(names(mRNA.train), affinities[[i]]), i] = 1
    }    
    model$fit.ridge = linearRidge(mRNA.train ~ ., data=data.frame(x.train))
  }
  else{
    warning("MAP configuration is the empty model!")   
    model$fit.ridge =lm(mRNA.train ~ 1)
  }  
  model$fit.ridge
}
