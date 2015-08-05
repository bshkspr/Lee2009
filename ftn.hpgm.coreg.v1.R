ftn.hpgm.coreg.v1 = function(spclst, chip.pval, hg.pval) {

##################################################################################

### !!! Inputs !!! ###
# 1. spclst    (list)   : a list of characters
# 2. chip.pval (number) : p-value threshold for ChIP-chip data (real value in [0,1])
# 3. hg.pval   (number) : p-value threshold for the hypergeometric score (real value in [0,1])

### !!! Outputs !!! ###
# 1. RESULTS (list) : a list of complexes with two objects, $pval and $allcrg.
#    (1) $pval   : a matrix giving significant co-regulating TFs and their HG p-value 
#                  at the given threshold (hg.pval)
#    (2) $allcrg : a list of significant co-regulating TFs and their target genes

### !!! Purpose !!! ###
# To identify co-regulated components for each TF in a complex using the hypergeometric 
# p-value in Simonis et al.'s paper (Simonis et al., Genome Biology 2004) 
# at a given ChIP-chip p-value threshold

### !!! Pitfalls !!! ###
# 1. This hypergeometric test is sensitive enough for small intersections, so that 
#    it is not appropriate if we want to test most of member proteins to be co-regulated.

### !!! Copyright !!! ###
# 1. Author   : H.F.L.
# 2. E-mail   : lee@molgen.mpg.de
# 3. Created  : 2006.06.09
# 4. Modified : 2009.02.10

#################################################################################

   os.win <- length(grep('w32', sessionInfo()$R.version$os)) != 0

#   if ( !exists('harb04.pval') ) {
#	if ( os.win ) load('./data/harb04.pval.RData')
#	else load('./data/harb04.pval.RData')
#   }

   time.start = date()

   all.compo = unique(unlist(spclst))
   all.balls = unique(c(rownames(harb04.pval), all.compo))
   sz.all.balls = length(all.balls)

   regulon = apply(harb04.pval, 2, function(x) 
	rownames(harb04.pval)[which(x < chip.pval)])

   regulon = regulon[!sapply(regulon, is.null)]

   RESULTS = lapply(spclst, as.null)

   for ( i in 1:length(spclst) ) {
	sltd.all = spclst[[i]]
	sz.sltd.all = length(sltd.all)

	if ( sz.sltd.all < 2 ) next

	tmp1 = NULL
	tmp.tfs = NULL
	tmp.crg = lapply(regulon, as.null)
	for ( j in 1:length(regulon) ) {
		wh.balls = regulon[[j]]
		sz.wh.balls = length(wh.balls)

		sltd.wh = sltd.all[sltd.all %in% wh.balls]
		sz.sltd.wh = length(sltd.wh)

		if ( sz.sltd.wh < 2 ) next

		tmp.crg[[j]] = sltd.wh

		# >>> Hypergeometric p-value (HPval) <<< #
		#     : P[X >= K] == 1 - P[X < K] == 1 - P[X <= (K-1)]
		hyp.pval = phyper(sz.sltd.wh - 1, sz.wh.balls, 
			sz.all.balls - sz.wh.balls, sz.sltd.all, lower.tail = F)

		tmp.summ = c(sz.wh.balls, sz.sltd.all, sz.sltd.wh, sz.all.balls, hyp.pval)
		tmp1 = rbind(tmp1, tmp.summ)

		tmp.tfs = c(tmp.tfs, names(regulon)[j])
	}

	if ( is.null(tmp1) ) next

	colnames(tmp1) = c('Regulon', 'Cpx', 'Coreg', 'BG', 'HPval')
	rownames(tmp1) = tmp.tfs

	if ( nrow(tmp1) > 1 )
		tmp1 = tmp1[order(tmp1[, 'HPval']), ]

	tmp2 = which(tmp1[, 'HPval'] < hg.pval)

	if ( length(tmp2) == 0 ) next
	if ( length(tmp2) == 1 ) {
		tmp3 = t(tmp1[tmp2, ])
		rownames(tmp3) = rownames(tmp1)[tmp2]
		tmp1 = tmp3
	} else tmp1 = tmp1[tmp2, ]

	RESULTS[[i]]$pval = tmp1

	tmp.crg = tmp.crg[rownames(tmp1)]
	RESULTS[[i]]$allcrg  = tmp.crg
   }

   RESULTS = RESULTS[sapply(RESULTS, length) != 0]

   time.end = date()
   cat('Program Start :', time.start, '\n')
   cat('Program End   :', time.end, '\n\n')
   cat('Done!\n\n')

   return(RESULTS)
}

