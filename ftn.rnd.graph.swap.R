ftn.rnd.graph.swap = function(Input1.graph, Input2.no.swaps) {

##################################################################################

### !!! Inputs !!! ###
# 1. Input1.graph (matrix) : a bipartite graph or network of proteins and complexes.
#    The format of the input graph (Input1.graph) should be N by 2, where N is the 
#    number of edges and the two columns correspond to two kinds of vertices (protein 
#    vertex and complex vertex) connected to each edge.
#    > For example, 
#         Input1.graph = as.matrix(read.table('c:/research/028_tomasz/data/gavin06.cmplx.pairs.txt'))
# 2. Input2.no.swaps (integer) : a number of random edge swaps
#    > For example, 
#         Input2.no.swaps = nrow(Input1.graph) * 1

### !!! Outputs !!! ###
# 1. rnd.graph (matrix) : a matrix of a randomized graph with the same format as 
#    the input graph (Input1.graph)

### !!! Purpose !!! ###
# To generate a simulated graph for a given graph (Input1.graph) by a given number of 
# random edge swaps (Input2.no.swaps)

### !!! Notes !!! ###
# 1. The randomization takes about 30 seconds for 1000 swaps.

### !!! Pitfalls !!! ###

### !!! Copyright !!! ###
# 1. Author   : H.F.L.
# 2. E-mail   : lee@molgen.mpg.de
# 3. Created  : 2006.06.17
# 4. Modified : 2006.06.20

#################################################################################

   start.time = date()
   runtime.st = proc.time()

   colnames(Input1.graph) = c('Prt', 'Cpx')
   rownames(Input1.graph) = paste('edge', 1:nrow(Input1.graph), sep='')

   rnd.graph = as.matrix(Input1.graph)
   # prt.set = unique(Input1.graph[, 'Prt'])
   # cpx.set = unique(Input1.graph[, 'Cpx'])
   edges.update = rownames(rnd.graph)

   for ( i in 1:Input2.no.swaps ) {
	# >>> 1. Select a random seed edge connected to a protein vertex and a complex vertex.
	#        This selection is done from previously unswapped edges, if any, to minimize the 
	#        possible bias of swapping for a specific edge.
	if ( length(edges.update) == 1 ) { seed.edge = edges.update
	} else seed.edge = sample(edges.update, 1)
	seed.prt  = rnd.graph[seed.edge, 'Prt']
	seed.cpx  = rnd.graph[seed.edge, 'Cpx']

	# >>> 2. Identify all complex/protein vertices to which the seed protein/complex vertex is linked.
	all.cpx.for.seed.prt = rnd.graph[rnd.graph[, 'Prt'] %in% seed.prt, 'Cpx']
	all.prt.for.seed.cpx = rnd.graph[rnd.graph[, 'Cpx'] %in% seed.cpx, 'Prt']

	# >>> 3. Identify all edges for which the seed protein/complex vertex is NOT relevant at all. 
	#        The identified edges are possbile edges among which one is to be swapped with the seed 
	#        edge. This is because the MAXIMUM number of possible edges between EACH protein vertex 
	#        and EACH complex vertex is ONE (i.e., no multi-edges).
	edges.pssb = rnd.graph
	edges.pssb = edges.pssb[!edges.pssb[, 'Prt'] %in% all.prt.for.seed.cpx, ]
	edges.pssb = edges.pssb[!edges.pssb[, 'Cpx'] %in% all.cpx.for.seed.prt, ]
	if ( !is.matrix(edges.pssb) ) {
		edges.pssb = t(edges.pssb)
		rownames(edges.pssb) = 
			names(which(apply(rnd.graph, 1, function(x) sum(x == edges.pssb) == 2)))
	}

	# sum(edges.pssb[, 'Prt'] %in% all.prt.for.seed.cpx)
	# [1] 0
	# sum(edges.pssb[, 'Cpx'] %in% all.cpx.for.seed.prt)
	# [1] 0

	# >>> 4. Select a random edge with which the seed edge is to be swapped.
	if ( length(rownames(edges.pssb)) == 1 ) { swap.edge = rownames(edges.pssb)
	} else swap.edge = sample(rownames(edges.pssb), 1)

	swap.cpx  = rnd.graph[swap.edge, 'Cpx']
	# swap.prt  = rnd.graph[swap.edge, 'Prt']

	# >>> 5. BEFORE SWAPPING, check and remember all protein vertices for each of the seed and 
	#        swap complex vertices to avoid the impermissible swapping.
	#        This is compared with AFTER SWAPPING below (step 7).
	before.seed = rnd.graph[rnd.graph[, 'Cpx'] %in% seed.cpx, 'Prt']
	before.swap = rnd.graph[rnd.graph[, 'Cpx'] %in% swap.cpx, 'Prt']

	# >>> 6. Swap the random two edges by exchanging the seed complex vertex with the random complex vertex
	rnd.graph[seed.edge, 'Cpx'] = swap.cpx
	rnd.graph[swap.edge, 'Cpx'] = seed.cpx

	# >>> 7. AFTER SWAPPING, check and remember all protein vertices for each of the seed and 
	#        swap complex vertices to avoid the impermissible swapping.
	#        This is compared with BEFORE SWAPPING above (step 5).
	after.seed = rnd.graph[rnd.graph[, 'Cpx'] %in% seed.cpx, 'Prt']
	after.swap = rnd.graph[rnd.graph[, 'Cpx'] %in% swap.cpx, 'Prt']

	# >>> 8. Check if the swapped complexes happen to be (topologically) the same as before.
	#        If so, undo the swapping.
	check.redund = setequal(before.seed, after.swap) & setequal(before.swap, after.seed)
	if ( check.redund ) {
		rnd.graph[seed.edge, 'Cpx'] = seed.cpx
		rnd.graph[swap.edge, 'Cpx'] = swap.cpx
	}

	# >>> 9. If the swapping has not been redundant, update the possbile edges 
	#        to select for the swapping at the next iteration.
	#        This is to ensure that all edges in the input graph are swapped at least once.
	if ( !check.redund ) edges.update = setdiff(edges.update, c(seed.edge, swap.edge))
	if ( length(edges.update) == 0 ) edges.update = rownames(rnd.graph)
   }

   # >>> 10. Check the degree distributions for protein and complex vertices
   deg.prt = table(Input1.graph[, 'Prt'])
   deg.cpx = table(Input1.graph[, 'Cpx'])

   deg.rnd.prt = table(rnd.graph[, 'Prt'])
   deg.rnd.cpx = table(rnd.graph[, 'Cpx'])

   deg.check = c(Prt = !identical(deg.prt, deg.rnd.prt), Cpx = !identical(deg.cpx, deg.rnd.cpx))
   if ( sum(deg.check) != 0 ) {
	if ( deg.check['Prt'] ) cat('The protein-degree distribution is not preserved!\n')
	if ( deg.check['Cpx'] ) cat('The complex-degree distribution is not preserved!\n')
	stop('The randomization of the input graph has failed!\n\n')
   }

   # >>> 11. Check identical complexes in the input graph (Input1.graph) and the randomized graph (rnd.graph)
   # ori.ntwk = split(as.character(Input1.graph[, 'Prt']), Input1.graph[, 'Cpx'])
   # rnd.ntwk = split(as.character(rnd.graph[, 'Prt']), rnd.graph[, 'Cpx'])

   # identi = sapply(ori.ntwk, function(x) sapply(rnd.ntwk, function(y) setequal(x, y)))
   # sum(identi)
   # which(identi, arr.ind=T)

   runtime.ed = proc.time()
   end.time = date()

   cat('Program Start :', start.time, '\n')
   cat('Program End   :', end.time, '\n')
   cat('Running Time  :', runtime.ed - runtime.st, '\n')
   cat('\nDone!\n\n')

   return(rnd.graph)

}

