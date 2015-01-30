dmixture = function(x, para, lambda){
	lambda[1]*dbeta(x, para[1], 1) + (1 - lambda[1])*dbeta(x, 1, para[2])
}

pmixture = function(x, para, lambda){
	lambda[1]*pbeta(x, para[1], 1) + (1 - lambda[1])*pbeta(x, 1, para[2])
}

qmix <- function(p,a,lambda,nbisect=20){			
	n <- length(p)
	mid <- rep(0,n)
	top <- rep(1,n)
	bot <- rep(0,n)
	gohigher <- rep(FALSE,n)
	for (j in 1:nbisect){
		mid <- (top+bot)/2
		gohigher <- (pmixture(mid,a,lambda)<p)		
		bot[gohigher] <- mid[gohigher]
		top[!gohigher] <- mid[!gohigher]
	}
	return(mid)
}

rmix <- function(n,a,lambda){
	u<-runif(n)
	return(qmix(u,a,lambda))
}

getComponent = function(x, i, para){
	if(i == 1)
		return(dbeta(x, para[1], 1))
	else
		return(dbeta(x, 1, para[2]))
}

loglik = function(para, X, lambda){
	para = exp(para) + 1		
	-sum(log(dmixture(X, para, lambda)))
}

bmixture.mle = function(X, para, lambda){
	res <- optim(log(para) - 1, fn=loglik, X=X,lambda=lambda, method="L-BFGS-B", upper=5, lower=-10)	
	list(para=exp(res$par) + 1, lambda=lambda, logLik=res$value+1e-50)
}

beta.EM = function(X, alpha=1.5, beta=5, lambda=0.5, tol=1e-4){
	converged <- FALSE		
	para <- c(alpha, beta)		
	ncomp <- 2
	lambda = c(lambda, 1-lambda)
	cat("start values: (alpha, beta) = ", para, "lambda = ",lambda,"\n")
	logLik <- 0
	Z <- matrix(0,nrow=length(X),ncol=ncomp)
	while(!converged){	
# 	E-step		
		for(i in 1:ncomp)
			Z[,i] <- lambda[i]*getComponent(X,i,para)/dmixture(X, para, lambda)
		lambdanew <- apply(Z,2,mean, na.rm=TRUE)
# 	M-step		
		paras <- bmixture.mle(X,para,lambdanew)	
		
		converged <- (abs(logLik/paras$logLik - 1) < tol)
		para <- paras$para
		lambda <- lambdanew
		logLik <- paras$logLik						
		cat("logLik = ", logLik, "(alpha, beta) = ", para, "lambda = ",lambda,"\n")		
	}
	print("converged!")	
	hist(X,probability=TRUE, main="histogram")
	x <- 1:100/100
	lines(x,dmixture(x,para,lambda),lwd=3)
	lines(x,lambda[1]*dbeta(x, para[1], 1),lwd=2,col="red")	
	lines(x,lambda[2]*dbeta(x, 1, para[2]),lwd=2,col="green")
	return(list(para=para,lambda=lambda,Z=Z,logLik=logLik))
}

suggestThreshold = function(prob, fpr=0.001){
  em = beta.EM(prob)
  tau1 = qbeta(1 - fpr, 1, em$para[2])
 # tau2 = optim(0.5, fn=function(x) 0.5*(log(beta(1, em$para[2])) - log(beta(em$para[1], 1)) + (em$para[1] - 1)*log(x) - (em$para[2] - 1)*log(1 - x) - lor)^2, 
  #      gr=function(x) 1/log(2)*(1/x * (em$para[1] - 1) - (em$para[2] - 1)/(x-1)), lower=0.01, upper=0.99)$par   
  #print(tau1)
  #print(tau2)
  #pmax(tau1, tau2)
 tau1
}
