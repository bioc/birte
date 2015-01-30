limmaAnalysis = function(dat, design=NULL, contrast) {
	if(class(dat) == "ExpressionSet")
		dat = exprs(dat)
	if(is.null(design)){
		design = model.matrix(~0 + factor(colnames(dat), levels=unique(colnames(dat))))
		colnames(design) <- unique(colnames(dat))		
	}		
	if(length(contrast) > 1)
		stop("limmaAnalysis supports only to extract one contrast!")
    fit = lmFit(dat, design)	
    contrast <- makeContrasts(contrasts=contrast, levels=design)
    fit2 <- contrasts.fit(fit, contrast)
    fit3 = eBayes(fit2)
    tab = topTable(fit3, number=Inf)	
    names(fit3$s2.post) = rownames(dat)
    list(pvalue.tab=tab, lm.fit=fit3, design=design)
}

