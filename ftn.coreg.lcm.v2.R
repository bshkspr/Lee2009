ftn.coreg.lcm.v2 = function(spclst, hgpval.regulon) {

##################################################################################

### !!! Inputs !!! ###
# 1. spclst (list) : a list of complexes or clusters
# 2. hgpval.regulon (list) : a list of matrices with hypergeometric p-values for each 
#                    regulon in each complex (identical format to the output of the 
#                    funtion 'ftn.hgpval.regulon.R')

### !!! Outputs !!! ###
# 1. OUTPUT (matrix) : local co-regulation measures, LCM, using the weighted sum of 
#                      the fractions of co-regulated components

### !!! Purpose !!! ###
# To calculate the local co-regulation measure (LCM) which is defined to be the 
# weighted sum of the fractions of co-regulated complex components by regulons 
# weighted by -log of the hypergeometric probabilities (p-value) of those 
# components at a given ChIP-chip p-value threshold in a population of complexes

### !!! Notes !!! ###
# (1) The weighted sum of LCM is converted into log score to scale down. 
#     The term 1 is added to that sum in order to make the log score >= 0.
# (2) This LCM definition yields the overlap of top scored complexes with 
#     functionally enriched complexes as good as the best overlap from LCM1 
#     ('ftn.coreg.lcm.v1') and better than the arithmetic mean of -log score 
#     at the hypergeometric p-value threshold 0.05.

### !!! Copyright !!! ###
# 1. Author   : H.F.L.
# 2. E-mail   : lee@molgen.mpg.de
# 3. Created  : 2006.06.30
# 4. Modified : 2006.08.27

#################################################################################

   # >>> To check if the two lists have identical names <<< #
   if ( !all(names(spclst) == names(hgpval.regulon)) ) 
	stop('Two lists do not seem to be from the same source! The order of the names must match.')

#   os.win <- length(grep('w32', sessionInfo()$R.version$os)) != 0

#   if ( !exists('') ) {
#	if ( os.win ) load('c:/research/028_tomasz/data/') else
#		load('/project/Scer/028_tomasz/data/')
#   }

#   time.start = date()
#   runtime.st = proc.time()

   # >>> DO NOT take into account for non-coregulated complexes in a population <<< #
   nocrg = which(sapply(hgpval.regulon, is.null))
   hgpval.regulon = hgpval.regulon[-nocrg]

   LCM.all = NULL

   for ( i in 1:length(hgpval.regulon) ) {
	tmp1 = hgpval.regulon[[i]]

	LCM.weights = tmp1[, 'Coreg'] / tmp1[, 'Cmplx']
	LCM.log.hgpval = -log(tmp1[, 'HGPval'])

	# >>> LCM, local co-regulation measure (LCM) for each complex 
	#     : log of the weighted sum of all co-regulated fractions added to 1 <<< #
	LCM.each = log(1 + sum(LCM.weights * LCM.log.hgpval))

	# >>> LCM2s and corresponding weights for all complexes <<< #
	LCM.all  = c(LCM.all, LCM.each)
   }

   names(LCM.all) = names(hgpval.regulon)

#   time.end = date()
#   cat('Program Start :', time.start, '\n')
#   cat('Program End   :', time.end, '\n\n')
#   cat('Running Time  :', (proc.time() - runtime.st)[3], 'sec\n\n')

   cat('Done!\n\n')

   return(LCM.all)
}

