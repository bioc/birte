getPotentialSwaps = function(genesets, perc.overlap.cutoff = 0.8, integer.id=TRUE) {
	T_potential_swaps = NULL
	if(! is.null(genesets) && length(genesets) > 0) {		
		tf = genesets
		mtf = matrix(nrow=length(tf), ncol=length(tf))		
		for(i in 1:length(tf)) {
			if(i < length(tf)){
				for(j in (i+1):length(tf)) {				
					mtf[i,j] = length(intersect(tf[[i]], tf[[j]]))/ length(union(tf[[i]], tf[[j]]))			
					mtf[j,i] = mtf[i,j]
				}
			}
		}

		sel_tfswaps = mtf > perc.overlap.cutoff
		T_potential_swaps = sapply(1:length(tf), function(x) {names(tf)[sel_tfswaps[x,]]})
		names(T_potential_swaps) = names(tf)
		if(integer.id) {
			T_potential_swaps = lapply(T_potential_swaps, function(x) {which(names(tf) %in% x)})
		}		
	}
	T_potential_swaps
}

