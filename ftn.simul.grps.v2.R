ftn.simul.grps.v2 = function(spcl) {

##################################################################################

### !!! Inputs !!! ###
# 1. spcl (list or matrix) : a group of clusters (e.g., complexes) in either 
# 			     a list or matrix format
# 	> An example list of clusters, 
# 	  : spcl = list(cx1 = c('p1', 'p2', 'X'), 
#			cx2 = c('p1', 'p3', 'X'), 
#			cx3 = c('p3', 'p4', 'p5', 'p6', 'X'))
# 	> An example matrix of clusters, 
# 	  : spcl = rbind(cbind('cx1', c('p1', 'p2', 'X')), 
#			 cbind('cx2', c('p1', 'p3', 'X')), 
#			 cbind('cx3', c('p3', 'p4', 'p5', 'p6', 'X')))

### !!! Outputs !!! ###
# 1. OUTPUT (list) : a list of random clusters by shuffling elements

### !!! Purpose !!! ###
# To simulate a given list of clusters (the input, 'spcl') by shuffling their 
# elements with multiplicities >= 1

### !!! Pitfalls !!! ###

### !!! Copyright !!! ###
# 1. Author   : H.F.L.
# 2. E-mail   : lee@molgen.mpg.de
# 3. Created  : 2006.07.05
# 4. Modified : 2006.07.05

#################################################################################

#   os.win <- length(grep('w32', sessionInfo()$R.version$os)) != 0

#   if ( !exists('......') ) {
#	if ( os.win ) load('c:/research/028_tomasz/......')
#	else load('/project/Scer/028_tomasz/......')
#   }

   time.start = date()
   runtime.st = proc.time()

   # >>> To name the input clusters for convenience <<< #
   # names(spcl) = paste('cls', 1:length(spcl), sep='')

   # >>> To convert the input matrix of clusters to a list format <<< #
   if ( is.matrix(spcl) ) spcl = split(spcl[, 2], spcl[, 1])

   # >>> Output : a randomized list of clusters (e.g. multiprotein complexes) <<< #
   OUTPUT = lapply(spcl, as.null)

   # >>> Unique elements and their multiplicities in the clusters (in the decreasing order) <<< #
   foo1 = sort(table(unlist(spcl)), decreasing = T)

   # >>> Sizes of the clusters <<< #
   foo2 = sapply(spcl, length)

   # >>> To distribute each element to random clusters according to its multiplicity (in the decreasing order) <<< #
   for ( i in 1:length(foo1) ) {

	# >>> To check all possible clusters for the distribution <<< #
	foo3 = names(which(foo2 > 0))

	# >>> To randomly select clusters as many as the multiplicity of the element under distribution <<< #
	if ( length(foo3) > 1 ) foo4 = sample(foo3, foo1[i]) else foo4 = foo3

	# >>> To avoid identical clusters <<< #
	repeat {
		z1 = names(which(foo2[foo4] == 1))
		z2 = which(duplicated(OUTPUT[z1]))
		if ( length(z2) == 0 ) break else {
			foo3 = setdiff(foo3, foo4)
			foo4 = setdiff(foo4, z1[z2])
			if ( length(foo3) > 1 ) 
				foo4 = c(foo4, sample(foo3, foo1[i]-length(foo4))) else foo4 = c(foo4, foo3)
		}
	}

	# >>> To distribute the element to the selected clusters <<< #
	for ( j in 1:length(foo4) ) {
		OUTPUT[[foo4[j]]] = c(OUTPUT[[foo4[j]]], names(foo1)[i])

		# >>> To update the sizes of all clusters for the next distribution of an element <<< #
		foo2[foo4[j]] = foo2[foo4[j]] - 1
	}
   }

   # names(OUTPUT) = paste(names(spcl), 's', sep='')

   # >>> To confirm the uniqueness of all the simulated clusters <<< #
   # chk.unq = sapply(OUTPUT, function(x) sapply(OUTPUT, function(y) setequal(x, y)))
   # diag(chk.unq) = FALSE
   # zzz[lower.tri(chk.unq)] = FALSE
   # if ( sum(chk.unq) != 0 ) stop('The randomization has generated identical clusters!\n')

   time.end = date()
   cat('Program Start :', time.start, '\n')
   cat('Program End   :', time.end, '\n\n')
   cat('Running Time  :', (proc.time() - runtime.st)[3], 'sec\n\n')
   cat('Done!\n\n')

   return(OUTPUT)
}

