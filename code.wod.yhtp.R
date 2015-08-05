
# >>> Shell command <<< #
# nohup echo "" | cat - code.wod.yhtp.R | nice R --vanilla --slave &

################################################################################

### !!! Title !!! ###
# "SYSTEMS-LEVEL TRANSCRIPTIONAL CO-REGULATION OF YEAST PROTEIN COMPLEXES"

### !!! Copyright !!! ###
# 1. Author   : J.W. Lee
# 2. E-mail   : lee55@fas.harvard.edu
# 3. Created  : 2006.09.16
# 4. Modified : 2009.02.13

################################################################################


# >>> PARAMETERS <<< #
DATA_FILE_1 = '/home/hojoon/projects/028_snr/data/wodak2008b_YORF.txt'
DATA_FILE_2 = '/home/hojoon/projects/028_snr/data/harb04.pval.RData'
FILENAME_0 = 'wod.yhtp'

RND_SIMUL_RANGE = 1:1000
ERR_SIMUL_RANGE = 1:1000

CHIP_THRS = c(0.001, 0.005, seq(0.01, 0.05, 0.005))
ERROR_RATES = seq(5, 20, 5)

DATASET_DIR = 'wod/'
CLUSTER_DIR = 'yhtp/'


### >>> FOR DIRECTORY LIST <<< ###
CURRENT_DIR = '/home/hojoon/projects/028_snr/'
DATA_DIR = 'data/'
CODE_DIR = 'code/'
RESULTS_DIR = 'results/'

SIMUL_DIR = 'simul/'
SIMUL_ERROR_DIR = 'rbst/'
SIMUL_HPGM_DIR = 'hpgm.PCMshff/'

CHIP_THRS_DIR = paste('p', sub('^0[.]', '', CHIP_THRS), '/', sep='')

ANALYSIS_DIR = 'rbst/'
NETWORK_TYPE_DIR = 'PCM/'

ERROR_DIR = paste('err', ERROR_RATES, '/', sep='')

GCM_DIR = 'GCM/'
LCM_DIR = 'LCM/'


### >>> FOR DIRECTORY AND FILE NAMES <<< ###
NAME_HPGM_REAL_POPUL_DIR = paste(CURRENT_DIR, RESULTS_DIR, DATASET_DIR, CLUSTER_DIR, GCM_DIR, sep='')
NAME_HPGM_REAL_POPUL_FILE_1 = paste(FILENAME_0, '.p', sep='')
NAME_HPGM_REAL_POPUL_FILE_2 = '.hpgm.RData'

NAME_LCM_REAL_RESULT_DIR = paste(CURRENT_DIR, RESULTS_DIR, DATASET_DIR, CLUSTER_DIR, GCM_DIR, LCM_DIR, sep='')
NAME_LCM_REAL_RESULT_REAL_FILE = paste(FILENAME_0, '.LCM2.diff.chip.thrs.real.RData', sep='')
NAME_LCM_REAL_RESULT_WILCOXON_TEST_FILE = paste(FILENAME_0, '.LCM2.diff.chip.thrs.mwtest.results.RData', sep='')
NAME_LCM_REAL_RESULT_WILCOXON_PVAL_FILE = paste(FILENAME_0, '.LCM2.diff.chip.thrs.mwtest.pval.txt', sep='')
NAME_LAMBDA_REAL_RESULT_ZSC_PVAL_FILE = paste(FILENAME_0, '.LAMBDA.arimean.LCM2.diff.chip.thrs.zsc.pval.txt', sep='')
NAME_LAMBDA_REAL_RESULT_ALL_FILE = paste(FILENAME_0, '.LAMBDA.arimean.LCM2.diff.chip.thrs.RData', sep='')

NAME_RND_POPUL_DIR  = paste(CURRENT_DIR, DATA_DIR, SIMUL_DIR, DATASET_DIR, CLUSTER_DIR, sep='')
NAME_RND_POPUL_FILE = paste(FILENAME_0, '.PCMshff.s', sep='')

NAME_HPGM_RND_POPUL_DIR = paste(CURRENT_DIR, RESULTS_DIR, DATASET_DIR, CLUSTER_DIR, SIMUL_DIR, SIMUL_HPGM_DIR, sep='')
NAME_HPGM_RND_POPUL_FILE_1 = paste(FILENAME_0, '.p', sep='')
NAME_HPGM_RND_POPUL_FILE_2 = '.hpgm.PCMshff.s'

NAME_LCM_RND_POPUL_DIR = paste(CURRENT_DIR, RESULTS_DIR, DATASET_DIR, CLUSTER_DIR, GCM_DIR, LCM_DIR, SIMUL_DIR, sep='')
NAME_LCM_RND_POPUL_FILE_0 = paste(FILENAME_0, '.chip_thrs_p', sep='')
NAME_LCM_RND_POPUL_FILE_1 = '.LCM2.PCMshff.simul1000.txt'
NAME_LCM_RND_POPUL_FILE_2 = '.LCM2.PCMshff.simul1000.RData'

NAME_ERR_POPUL_DIR  = paste(CURRENT_DIR, DATA_DIR, SIMUL_DIR, SIMUL_ERROR_DIR, NETWORK_TYPE_DIR, ERROR_DIR, DATASET_DIR, CLUSTER_DIR, sep='')
NAME_ERR_POPUL_FILE = paste(FILENAME_0, '.err', ERROR_RATES, '.PCMshff.s', sep='')

NAME_HPGM_ERR_POPUL_DIR = paste(CURRENT_DIR, RESULTS_DIR, ANALYSIS_DIR, NETWORK_TYPE_DIR, ERROR_DIR, DATASET_DIR, CLUSTER_DIR, SIMUL_DIR, SIMUL_HPGM_DIR, sep='')
NAME_HPGM_ERR_POPUL_FILE_1 = paste(FILENAME_0, '.err', ERROR_RATES, '.p', sep='')
NAME_HPGM_ERR_POPUL_FILE_2 = '.hpgm.PCMshff.s'

NAME_LCM_ERR_POPUL_DIR = paste(CURRENT_DIR, RESULTS_DIR, ANALYSIS_DIR, NETWORK_TYPE_DIR, ERROR_DIR, DATASET_DIR, CLUSTER_DIR, GCM_DIR, LCM_DIR, SIMUL_DIR, sep='')
NAME_LCM_ERR_POPUL_FILE_0 = paste(FILENAME_0, '.err', ERROR_RATES, '.sim1000.chip_thrs_p', sep='')
NAME_LCM_ERR_POPUL_FILE_1 = '.LCM2.txt'
NAME_LCM_ERR_POPUL_FILE_2 = '.LCM2.RData'

NAME_LCM_ERR_RESULT_DIR = paste(CURRENT_DIR, RESULTS_DIR, ANALYSIS_DIR, NETWORK_TYPE_DIR, ERROR_DIR, DATASET_DIR, CLUSTER_DIR, GCM_DIR, LCM_DIR, sep='')
NAME_LCM_ERR_RESULT_WILCOXON_TEST_FILE = paste(FILENAME_0, '.err', ERROR_RATES, '_sim1000.LCM2.diff.chip.thrs.PCMshff_sim1000.mwtest.results.RData', sep='')
NAME_LCM_ERR_RESULT_WILCOXON_PVAL_FILE = paste(FILENAME_0, '.err', ERROR_RATES, '_sim1000.LCM2.diff.chip.thrs.PCMshff_sim1000.mwtest.pval.txt', sep='')
NAME_LAMBDA_ERR_RESULT_ZSC_PVAL_FILE = paste(FILENAME_0, '.err', ERROR_RATES, '.LAMBDA.arimean.LCM2.diff.chip.thrs.zsc.pval.txt', sep='')
NAME_LAMBDA_ERR_RESULT_ALL_FILE = paste(FILENAME_0, '.err', ERROR_RATES, '.LAMBDA.arimean.LCM2.diff.chip.thrs.RData', sep='')


### >>> FOR DIRECTORY CREATION AND MANAGEMENT <<< ###
setwd(CURRENT_DIR)
if ( !RESULTS_DIR %in% paste(dir(), '/', sep='') )
	dir.create(RESULTS_DIR)
if ( !DATA_DIR %in% paste(dir(), '/', sep='') )
	dir.create(DATA_DIR)

setwd(paste(CURRENT_DIR, DATA_DIR, sep=''))
if ( !SIMUL_DIR %in% paste(dir(), '/', sep='') )
	dir.create(SIMUL_DIR)

setwd(paste(CURRENT_DIR, DATA_DIR, SIMUL_DIR, sep=''))
if ( !DATASET_DIR %in% paste(dir(), '/', sep='') )
	dir.create(DATASET_DIR)
if ( !SIMUL_ERROR_DIR %in% paste(dir(), '/', sep='') )
	dir.create(SIMUL_ERROR_DIR)

setwd(paste(CURRENT_DIR, DATA_DIR, SIMUL_DIR, DATASET_DIR, sep=''))
if ( !CLUSTER_DIR %in% paste(dir(), '/', sep='') )
	dir.create(CLUSTER_DIR)

setwd(paste(CURRENT_DIR, DATA_DIR, SIMUL_DIR, SIMUL_ERROR_DIR, sep=''))
if ( !NETWORK_TYPE_DIR %in% paste(dir(), '/', sep='') )
	dir.create(NETWORK_TYPE_DIR)

for ( i in 1:length(ERROR_RATES) ) {
   setwd(paste(CURRENT_DIR, DATA_DIR, SIMUL_DIR, SIMUL_ERROR_DIR, NETWORK_TYPE_DIR, sep=''))
   if ( !ERROR_DIR[i] %in% paste(dir(), '/', sep='') )
	dir.create(ERROR_DIR[i])

   setwd(ERROR_DIR[i])
   if ( !DATASET_DIR %in% paste(dir(), '/', sep='') )
	dir.create(DATASET_DIR)

   setwd(DATASET_DIR)
   if ( !CLUSTER_DIR %in% paste(dir(), '/', sep='') )
	dir.create(CLUSTER_DIR)
}
rm(i)

setwd(paste(CURRENT_DIR, RESULTS_DIR, sep=''))
if ( !ANALYSIS_DIR %in% paste(dir(), '/', sep='') )
	dir.create(ANALYSIS_DIR)

setwd(paste(CURRENT_DIR, RESULTS_DIR, ANALYSIS_DIR, sep=''))
if ( !NETWORK_TYPE_DIR %in% paste(dir(), '/', sep='') )
	dir.create(NETWORK_TYPE_DIR)

for ( i in 1:length(ERROR_RATES) ) {
   setwd(paste(CURRENT_DIR, RESULTS_DIR, ANALYSIS_DIR, NETWORK_TYPE_DIR, sep=''))
   if ( !ERROR_DIR[i] %in% paste(dir(), '/', sep='') )
	dir.create(ERROR_DIR[i])

   setwd(ERROR_DIR[i])
   if ( !DATASET_DIR %in% paste(dir(), '/', sep='') )
	dir.create(DATASET_DIR)

   setwd(DATASET_DIR)
   if ( !CLUSTER_DIR %in% paste(dir(), '/', sep='') )
	dir.create(CLUSTER_DIR)

   setwd(CLUSTER_DIR)
   if ( !SIMUL_DIR %in% paste(dir(), '/', sep='') )
	dir.create(SIMUL_DIR)

   setwd(SIMUL_DIR)
   if ( !SIMUL_HPGM_DIR %in% paste(dir(), '/', sep='') )
	dir.create(SIMUL_HPGM_DIR)

   setwd('../')
   if ( !GCM_DIR %in% paste(dir(), '/', sep='') )
	dir.create(GCM_DIR)

   setwd(GCM_DIR)
   if ( !LCM_DIR %in% paste(dir(), '/', sep='') )
	dir.create(LCM_DIR)

   setwd(LCM_DIR)
   if ( !SIMUL_DIR %in% paste(dir(), '/', sep='') )
	dir.create(SIMUL_DIR)
}
rm(i)

setwd(paste(CURRENT_DIR, RESULTS_DIR, sep=''))
if ( !DATASET_DIR %in% paste(dir(), '/', sep='') )
	dir.create(DATASET_DIR)

setwd(paste(CURRENT_DIR, RESULTS_DIR, DATASET_DIR, sep=''))
if ( !CLUSTER_DIR %in% paste(dir(), '/', sep='') )
	dir.create(CLUSTER_DIR)

setwd(paste(CURRENT_DIR, RESULTS_DIR, DATASET_DIR, CLUSTER_DIR, sep=''))
if ( !GCM_DIR %in% paste(dir(), '/', sep='') )
	dir.create(GCM_DIR)

setwd(paste(CURRENT_DIR, RESULTS_DIR, DATASET_DIR, CLUSTER_DIR, GCM_DIR, sep=''))
if ( !LCM_DIR %in% paste(dir(), '/', sep='') )
	dir.create(LCM_DIR)

setwd(paste(CURRENT_DIR, RESULTS_DIR, DATASET_DIR, CLUSTER_DIR, GCM_DIR, LCM_DIR, sep=''))
if ( !SIMUL_DIR %in% paste(dir(), '/', sep='') )
	dir.create(SIMUL_DIR)

setwd(paste(CURRENT_DIR, RESULTS_DIR, DATASET_DIR, CLUSTER_DIR, sep=''))
if ( !SIMUL_DIR %in% paste(dir(), '/', sep='') )
	dir.create(SIMUL_DIR)

setwd(paste(CURRENT_DIR, RESULTS_DIR, DATASET_DIR, CLUSTER_DIR, SIMUL_DIR, sep=''))
if ( !SIMUL_HPGM_DIR %in% paste(dir(), '/', sep='') )
	dir.create(SIMUL_HPGM_DIR)

setwd(paste(CURRENT_DIR, RESULTS_DIR, DATASET_DIR, CLUSTER_DIR, SIMUL_DIR, SIMUL_HPGM_DIR, sep=''))
for ( i in 1:length(CHIP_THRS_DIR) ) {
   if ( !CHIP_THRS_DIR[i] %in% paste(dir(), '/', sep='') )
	dir.create(CHIP_THRS_DIR[i])
}

setwd(CURRENT_DIR)


### >>> REQUIRED FUNCTIONS <<< ###
NAME_FTN_DIR = paste(CURRENT_DIR, CODE_DIR, sep='')
FUNCTION_1 = 'ftn.hpgm.coreg.v1.R'
FUNCTION_2 = 'ftn.coreg.lcm.v2.R'
FUNCTION_3 = 'ftn.simul.grps.v2.R'
FUNCTION_4 = 'ftn.rnd.graph.swap.R'


################################################################################
# ============================================================================ #
# 				>>> PART 0 <<<				       #

# General settings
# ============================================================================ #

# >>> Required data <<< #
GRAPH = as.matrix(read.table(DATA_FILE_1))
GRAPH_LIST = split(as.character(GRAPH[, 2]), GRAPH[, 1])

# >>> Required functions <<< #
source(paste(NAME_FTN_DIR, FUNCTION_1, sep=''))
source(paste(NAME_FTN_DIR, FUNCTION_2, sep=''))
source(paste(NAME_FTN_DIR, FUNCTION_3, sep=''))
source(paste(NAME_FTN_DIR, FUNCTION_4, sep=''))

# >>> GLOBAL VARIABLES <<< #
COLNAMES_REAL = c('ChIP_Pval', 'Real', 'Mean_Rnd', 'SD_Rnd', 'Z_score', 'P_value')
COLNAMES_ERR  = c('ChIP_Pval', 'Real_Err', 'Mean_Rnd', 'SD_Rnd', 'Z_score', 'P_value')

COMMAND_TAR_RND_POPUL = paste('tar zcvf ', FILENAME_0, '.sim1000.PCMshff.txt.tar *PCMshff.s*.txt', sep='')
COMMAND_TAR_RND_CHIP = paste('tar zcvf ', FILENAME_0, '.p', sep='')
COMMAND_TAR_HPGM = '.hpgm.sim1000.PCMshff.RData.tar *PCMshff.s*.RData'
COMMAND_TAR_ERR_0 = paste('tar zcvf ', FILENAME_0, '.err', sep='')
COMMAND_TAR_ERR_POPUL = '.sim1000.PCMshff.txt.tar *PCMshff.s*.txt'

start.time = date()


################################################################################
# ============================================================================ #
# 				>>> PART I <<<				       #

# To calculate the hypergeometric p-values for all regulons in each complex 
# and identify all co-regulated protein components 
# (at each ChIP-chip p-value threshold)
# : FUNCTION_1
# ============================================================================ #

setwd(NAME_HPGM_REAL_POPUL_DIR)

load(DATA_FILE_2)

for ( thrs in CHIP_THRS ) {
   cat('> PART 1/12 ; Threshold', thrs, '\n')

   tmp.thrs = sub('^0[.]', '', thrs)

   hpgm.real = ftn.hpgm.coreg.v1(GRAPH_LIST, thrs, 1)

   save(hpgm.real, file=paste(NAME_HPGM_REAL_POPUL_FILE_1, tmp.thrs, 
	NAME_HPGM_REAL_POPUL_FILE_2, sep=''))

}
rm(thrs)

cat('... End of PART I ...\n\n')
setwd(CURRENT_DIR)


################################################################################
# ============================================================================ #
# 				>>> PART II <<<				       #

# To calculate the local co-regulation measures, LCM_2, for the population 
# of real complexes
# : FUNCTION_2 (weighted sum of the fractions of co-regulated components)
# ============================================================================ #

setwd(NAME_LCM_REAL_RESULT_DIR)

lcm.real = vector('list', length = length(CHIP_THRS))
names(lcm.real) = CHIP_THRS

for ( thrs in CHIP_THRS ) {
   cat('> PART 2/12 ; Threshold', thrs, '\n')

   tmp.thrs = sub('^0[.]', '', thrs)

   load(paste(NAME_HPGM_REAL_POPUL_DIR, NAME_HPGM_REAL_POPUL_FILE_1, 
	tmp.thrs, NAME_HPGM_REAL_POPUL_FILE_2, sep=''))
   # 'hpgm.real'

   if ( length(hpgm.real) == 0 ) tmpsim.out = 0 else {
	for ( k in 1:length(hpgm.real) ) {
		colnames(hpgm.real[[k]]$pval)[c(2, 5)] = c('Cmplx', 'HGPval')
	}
	rm(k)

	hgpval.real = lapply(GRAPH_LIST, as.null)
	hgpval.real[names(hpgm.real)] = lapply(hpgm.real, function(x) x$pval)

	tmp.out = ftn.coreg.lcm.v2(GRAPH_LIST, hgpval.real)
   }

   lcm.real[[as.character(thrs)]] = sort(tmp.out, decreasing = T)

}
rm(thrs)

save(lcm.real, file=NAME_LCM_REAL_RESULT_REAL_FILE)

cat('... End of PART II ...\n\n')
setwd(CURRENT_DIR)


################################################################################
# ============================================================================ #
# 				>>> PART III <<<			       #

# To generate populations of randomized complexes by PCM shuffling
# : FUNCTION_3
# ============================================================================ #

setwd(NAME_RND_POPUL_DIR)

for ( i in RND_SIMUL_RANGE ) {
	cat('> PART 3/12 ; Simulation', i, '\n')

	tmpsim = ftn.simul.grps.v2(GRAPH)
	names(tmpsim) = paste(names(tmpsim), 's', sep='')

	pairs = NULL
	for ( j in names(tmpsim) ) {
		pairs = rbind(pairs, cbind(j, tmpsim[[j]]))
	}
	rm(j)

	write.table(pairs, file = paste(NAME_RND_POPUL_FILE, i, '.txt', sep=''), 
		quote=F, row.names=F, col.names=F, sep='\t')
}
rm(i)

system(COMMAND_TAR_RND_POPUL)

cat('... End of PART III ...\n\n')
setwd(CURRENT_DIR)


################################################################################
# ============================================================================ #
# 				>>> PART IV <<<				       #

# To calculate the hypergeometric p-value for each regulon in each randomized 
# complex and to identify all co-regulated protein components
# : FUNCTION_1
# ============================================================================ #

setwd(NAME_HPGM_RND_POPUL_DIR)

load(DATA_FILE_2)

for ( thrs in CHIP_THRS ) {

   tmp.thrs = sub('^0[.]', '', thrs)

   if ( !paste('p', tmp.thrs, sep='') %in% dir() )
	dir.create(paste('p', tmp.thrs, sep=''))

   setwd(paste('p', tmp.thrs, sep=''))

   for ( i in RND_SIMUL_RANGE ) {
	cat('> PART 4/12 ; Threshold', thrs, '; Simulation', i, '\n')

	FILENAME_PCF = paste(NAME_RND_POPUL_DIR, NAME_RND_POPUL_FILE, i, '.txt', sep='')
	tmpsim = read.table(FILENAME_PCF)
	tmpsim = split(as.character(tmpsim[, 2]), tmpsim[, 1])

	hpgm.sim = ftn.hpgm.coreg.v1(tmpsim, thrs, 1)

	save(hpgm.sim, file=paste(NAME_HPGM_RND_POPUL_FILE_1, tmp.thrs, 
		NAME_HPGM_RND_POPUL_FILE_2, i, '.RData', sep=''))
   }
   rm(i)

   system(paste(COMMAND_TAR_RND_CHIP, tmp.thrs, COMMAND_TAR_HPGM, sep=''))

   setwd('../')
}
rm(thrs)

cat('... End of PART IV ...\n\n')
setwd(CURRENT_DIR)


################################################################################
# ============================================================================ #
# 				>>> PART V <<<				       #

# To calculate the local co-regulation measure, LCM_2, for populations of 
# randomized complexes
# : FUNCTION_2 (weighted sum of the fractions of co-regulated components)
# ============================================================================ #

setwd(NAME_LCM_RND_POPUL_DIR)
FILENAME_PCF = paste(NAME_RND_POPUL_DIR, NAME_RND_POPUL_FILE, sep='')

for ( thrs in CHIP_THRS ) {
      tmp.thrs = sub('^0[.]', '', thrs)
      lcm.rnd = vector('list', length = RND_SIMUL_RANGE)

      for ( i in RND_SIMUL_RANGE ) {
	cat('> PART 5/12 ; Threshold', thrs, '; Simulation', i, '\n')

	tmp.file = paste(FILENAME_PCF, i, '.txt', sep='')
	tmp.sim = read.table(tmp.file)

	spcl.sim = split(as.character(tmp.sim[, 2]), tmp.sim[, 1])
# !!!	spcl.sim = split(as.character(tmp.sim[, 1]), tmp.sim[, 2])

	load(paste(NAME_HPGM_RND_POPUL_DIR, paste('p', tmp.thrs, '/', sep=''), 
		NAME_HPGM_RND_POPUL_FILE_1, tmp.thrs, NAME_HPGM_RND_POPUL_FILE_2, i, '.RData', sep=''))
	# "hpgm.sim"

	if ( length(hpgm.sim) == 0 ) tmpsim.out = 0 else {
	   for ( j in 1:length(hpgm.sim) ) {
		colnames(hpgm.sim[[j]]$pval)[c(2, 5)] = c('Cmplx', 'HGPval')
	   }
	   rm(j)

	   hgpval.sim = lapply(spcl.sim, as.null)
	   hgpval.sim[names(hpgm.sim)] = lapply(hpgm.sim, function(x) x$pval)

	   tmpsim.out = ftn.coreg.lcm.v2(spcl.sim, hgpval.sim)
	}

	cat(file=paste(NAME_LCM_RND_POPUL_FILE_0, tmp.thrs, NAME_LCM_RND_POPUL_FILE_1, sep=''), append=T, 
		tmpsim.out, sep='\n')

	lcm.rnd[[i]] = tmpsim.out
      }
      rm(i)

      save(lcm.rnd, file=paste(NAME_LCM_RND_POPUL_FILE_0, tmp.thrs, NAME_LCM_RND_POPUL_FILE_2, sep=''))
}
rm(thrs)

cat('... End of PART V ...\n\n')
setwd(CURRENT_DIR)


################################################################################
# ============================================================================ #
# 				>>> PART VI <<<				       #

# To calculate the Wilcoxon P-values for the population of real complexes with 
# respect to an ensemble of random populations
# ============================================================================ #

setwd(NAME_LCM_REAL_RESULT_DIR)

load(file = NAME_LCM_REAL_RESULT_REAL_FILE)
# "lcm.real"

mwtest = vector('list', length = length(CHIP_THRS))
names(mwtest) = CHIP_THRS

for ( thrs in CHIP_THRS ) {
	cat('> PART 6/12 ; Threshold', thrs, '\n')

	# >>> For the real data <<< #
	foo.real = lcm.real[[as.character(thrs)]]

	# >>> For the random data <<< #
	foo.rnd = scan(paste(NAME_LCM_RND_POPUL_DIR, NAME_LCM_RND_POPUL_FILE_0, 
		sub('^0[.]', '', thrs), NAME_LCM_RND_POPUL_FILE_1, sep=''))

	# >>> Two sample Wilcoxon rank sum test (Mann-Whitney) without continuity correction <<< #
	tmp.mwtest = wilcox.test(foo.real, foo.rnd, alternative='g', correct=F)
	mwtest[[as.character(thrs)]] = tmp.mwtest
}
rm(thrs)

mwtest.pval = as.matrix(sapply(mwtest, function(tmp) tmp$p.value))
rownames(mwtest.pval) = paste('chip_', rownames(mwtest.pval), sep='')

save(mwtest, file = NAME_LCM_REAL_RESULT_WILCOXON_TEST_FILE)
write.table(mwtest.pval, file = NAME_LCM_REAL_RESULT_WILCOXON_PVAL_FILE, quote = F, col.names = F, sep='\t')

cat('... End of PART VI ...\n\n')
setwd(CURRENT_DIR)


################################################################################
# ============================================================================ #
# 				>>> PART VII <<<			       #

# To calculate the global co-regulation measure, LAMBDA, its Z-scores and 
# P-values for the population of real complexes
# ============================================================================ #

setwd(NAME_LCM_REAL_RESULT_DIR)

load(file = NAME_LCM_REAL_RESULT_REAL_FILE)
# "lcm.real"

Lambda = vector('list', length = length(CHIP_THRS))
names(Lambda) = CHIP_THRS

for ( thrs in CHIP_THRS ) {
   cat('> PART 7/12 ; Threshold', thrs, '\n')

   # >>> For the real data <<< #
   foo.real = lcm.real[[as.character(thrs)]]
   lambda.real = mean(foo.real)

   Lambda[[as.character(thrs)]]$real = lambda.real

   # >>> For the random data <<< #
   load(paste(NAME_LCM_RND_POPUL_DIR, NAME_LCM_RND_POPUL_FILE_0, 
	sub('^0[.]', '', thrs), NAME_LCM_RND_POPUL_FILE_2, sep=''))
   # [1] "lcm.rnd"

   lambda.rnd = sapply(lcm.rnd, mean)

   Lambda[[as.character(thrs)]]$rnd = lambda.rnd

   # >>> For Z-scores and P-values <<< #
   rnd.mean = mean(lambda.rnd)
   rnd.sd   = sd(lambda.rnd)

   z.scores = (lambda.real - rnd.mean) / rnd.sd
   p.values = sum(lambda.rnd >= lambda.real) / length(lambda.rnd)

   results = cbind(Real = lambda.real, Mean_Rnd = rnd.mean, SD_Rnd = rnd.sd, 
	Z_score = z.scores, P_value = p.values)

   Lambda[[as.character(thrs)]]$res = results
}

foo.lam = t(sapply(Lambda, function(x) round(x$res, 3)))
foo.lam = cbind(rownames(foo.lam), foo.lam)
colnames(foo.lam) = COLNAMES_REAL

write.table(foo.lam, file = NAME_LAMBDA_REAL_RESULT_ZSC_PVAL_FILE, quote=F, row.names=F, sep='\t')
save(Lambda, file = NAME_LAMBDA_REAL_RESULT_ALL_FILE)

cat('... End of PART VII ...\n\n')
setwd(CURRENT_DIR)


################################################################################
# ============================================================================ #
# 				>>> PART VIII <<<			       #

# To generate error-containing complexes in a population of real complexes by 
# randomly swapping the membership of protein components 
# (or edges in the PCF graph)
# : FUNCTION_4
# ============================================================================ #

for ( i in 1:length(ERROR_RATES) ) {

   setwd(NAME_ERR_POPUL_DIR[i])

   # >>> Parameters <<< #
   ERRORS = ERROR_RATES[i] / 100
   N_SWAPS = nrow(GRAPH) * ERRORS

   # >>> MAIN <<< #
   for ( j in ERR_SIMUL_RANGE ) {
	cat('> PART 8/12 ; Error Rate', ERROR_RATES[i], '; Simulation', j, '\n')
# !!!	tmpsim = ftn.rnd.graph.swap(GRAPH, N_SWAPS)
	tmpsim = ftn.rnd.graph.swap(GRAPH[, c(2,1)], N_SWAPS)

	tmpsim = tmpsim[, c(2,1)]
	write.table(tmpsim, file=paste(NAME_ERR_POPUL_FILE[i], j, '.txt', sep=''), 
		quote=F, row.names=F, col.names=F, sep='\t')
   }
   rm(j)

   system(paste(COMMAND_TAR_ERR_0, ERROR_RATES[i], COMMAND_TAR_ERR_POPUL, sep=''))

}
rm(i)

cat('... End of PART VIII ...\n\n')
setwd(CURRENT_DIR)


################################################################################
# ============================================================================ #
# 				>>> PART IX <<<				       #

# To examine populations of error-containing complexes for the calculation of 
# the hypergeometric p-value for each regulon and to identify all co-regulated 
# protein components
# : FUNCTION_1
# ============================================================================ #

load(DATA_FILE_2)

for ( i in 1:length(ERROR_RATES) ) {

   setwd(NAME_HPGM_ERR_POPUL_DIR[i])
   FILENAME_PCF = paste(NAME_ERR_POPUL_DIR[i], NAME_ERR_POPUL_FILE[i], sep='')

   for ( thrs in CHIP_THRS ) {

      tmp.thrs = sub('^0[.]', '', thrs)

      if ( !paste('p', tmp.thrs, sep='') %in% dir() )
	dir.create(paste('p', tmp.thrs, sep=''))

      setwd(paste('p', tmp.thrs, sep=''))

      for ( j in ERR_SIMUL_RANGE ) {
	cat('> PART 9/12 ; Error Rate', ERROR_RATES[i], '% ; Threshold', thrs, '; Simulation', j, '\n')

	tmpsim = read.table(paste(FILENAME_PCF, j, '.txt', sep=''))
	tmpsim = split(as.character(tmpsim[, 2]), tmpsim[, 1])
# !!!	tmpsim = split(as.character(tmpsim[, 1]), tmpsim[, 2])

	hpgm.sim = ftn.hpgm.coreg.v1(tmpsim, thrs, 1)

	save(hpgm.sim, file=paste(NAME_HPGM_ERR_POPUL_FILE_1[i], 
		tmp.thrs, NAME_HPGM_ERR_POPUL_FILE_2, j, '.RData', sep=''))
      }
      rm(j)

      system(paste(COMMAND_TAR_ERR_0, ERROR_RATES[i], '.p', tmp.thrs, COMMAND_TAR_HPGM, sep=''))

      setwd('../')
   }
   rm(thrs)
}
rm(i)

cat('... End of PART IX ...\n\n')
setwd(CURRENT_DIR)


################################################################################
# ============================================================================ #
# 				>>> PART X <<<				       #

# To calculate the local co-regulation measure, LCM_2, for populations of 
# error-containing complexes
# : FUNCTION_2 (weighted sum of the fractions of co-regulated components)
# ============================================================================ #

for ( i in 1:length(ERROR_RATES) ) {

   setwd(NAME_LCM_ERR_POPUL_DIR[i])
   FILENAME_PCF = paste(NAME_ERR_POPUL_DIR[i], NAME_ERR_POPUL_FILE[i], sep='')

   for ( thrs in CHIP_THRS ) {
      tmp.thrs = sub('^0[.]', '', thrs)
      lcm.err = vector('list', length = ERR_SIMUL_RANGE)

      for ( j in ERR_SIMUL_RANGE ) {
	cat('> PART 10/12 ; Error Rate', ERROR_RATES[i], '% ; Threshold', thrs, '; Simulation', j, '\n')

	tmp.file = paste(FILENAME_PCF, j, '.txt', sep='')
	tmp.sim = read.table(tmp.file)

	spcl.sim = split(as.character(tmp.sim[, 2]), tmp.sim[, 1])
# !!!	spcl.sim = split(as.character(tmp.sim[, 1]), tmp.sim[, 2])

	load(paste(NAME_HPGM_ERR_POPUL_DIR[i], paste('p', tmp.thrs, '/', sep=''), 
		NAME_HPGM_ERR_POPUL_FILE_1[i], tmp.thrs, NAME_HPGM_ERR_POPUL_FILE_2, j, '.RData', sep=''))
	# "hpgm.sim"

	if ( length(hpgm.sim) == 0 ) tmpsim.out = 0 else {
	   for ( k in 1:length(hpgm.sim) ) {
		colnames(hpgm.sim[[k]]$pval)[c(2, 5)] = c('Cmplx', 'HGPval')
	   }
	   rm(k)

	   hgpval.sim = lapply(spcl.sim, as.null)
	   hgpval.sim[names(hpgm.sim)] = lapply(hpgm.sim, function(x) x$pval)

	   tmpsim.out = ftn.coreg.lcm.v2(spcl.sim, hgpval.sim)
	}

	cat(file=paste(NAME_LCM_ERR_POPUL_FILE_0[i], tmp.thrs, NAME_LCM_ERR_POPUL_FILE_1, sep=''), append=T, 
		tmpsim.out, sep='\n')

	lcm.err[[j]] = sort(tmpsim.out, decreasing = T)
      }
      rm(j)

      save(lcm.err, file=paste(NAME_LCM_ERR_POPUL_FILE_0[i], tmp.thrs, NAME_LCM_ERR_POPUL_FILE_2, sep=''))
   }
   rm(thrs)
}
rm(i)

cat('... End of PART X ...\n\n')
setwd(CURRENT_DIR)


################################################################################
# ============================================================================ #
# 				>>> PART XI <<<				       #

# To calculate the Wilcoxon P-values for populations of error-containing 
# complexes with respect to an ensemble of random populations
# ============================================================================ #

for ( i in 1:length(ERROR_RATES) ) {

   setwd(NAME_LCM_ERR_RESULT_DIR[i])

   mwtest = vector('list', length = length(CHIP_THRS))
   names(mwtest) = CHIP_THRS

   for ( thrs in CHIP_THRS ) {
      # >>> For the random data <<< #
      foo.rnd = scan(paste(NAME_LCM_RND_POPUL_DIR, NAME_LCM_RND_POPUL_FILE_0, 
	sub('^0[.]', '', thrs), NAME_LCM_RND_POPUL_FILE_1, sep=''))

      # >>> For the error-containing real data <<< #
      load(file = paste(SIMUL_DIR, NAME_LCM_ERR_POPUL_FILE_0[i], 
	sub('^0[.]', '', thrs), NAME_LCM_ERR_POPUL_FILE_2, sep=''))
      # "lcm.err"

      for ( k in 1:length(lcm.err) ) {
	cat('> PART 11/12 ; Error Rate', ERROR_RATES[i], '% ; Threshold', thrs, '; Simulation', k, '\n')

	foo.real = lcm.err[[k]]

	# >>> Two sample Wilcoxon rank sum test (Mann-Whitney) without continuity correction <<< #
	tmp.mwtest = wilcox.test(foo.real, foo.rnd, alternative='g', correct=F)
	mwtest[[as.character(thrs)]][[k]] = tmp.mwtest
      }
      rm(k)
   }
   rm(thrs)

   mwtest.pval = sapply(mwtest, function(tmp) sapply(tmp, function(x) x$p.value))
   colnames(mwtest.pval) = paste('chip_', colnames(mwtest.pval), sep='')

   save(mwtest, file = NAME_LCM_ERR_RESULT_WILCOXON_TEST_FILE[i])
   write.table(mwtest.pval, file = NAME_LCM_ERR_RESULT_WILCOXON_PVAL_FILE[i], quote = F, row.names = F, sep='\t')
}

cat('... End of PART XI ...\n\n')
setwd(CURRENT_DIR)


################################################################################
# ============================================================================ #
# 				>>> PART XII <<<			       #

# To calculate the global co-regulation measure, LAMBDA, its Z-scores and 
# P-values for populations of error-containing complexes
# ============================================================================ #

for ( i in 1:length(ERROR_RATES) ) {

   setwd(NAME_LCM_ERR_RESULT_DIR[i])

   Lambda = vector('list', length = length(CHIP_THRS))
   names(Lambda) = CHIP_THRS

   for ( thrs in CHIP_THRS ) {

      # >>> For the error-containing real data <<< #
      lambda.real = NULL
      for ( j in ERR_SIMUL_RANGE ) {
	cat('> PART 12/12 ; Error Rate', ERROR_RATES[i], '% ; Threshold', thrs, '; Simulation', j, '\n')

	load(file = paste(SIMUL_DIR, NAME_LCM_ERR_POPUL_FILE_0[i], 
		sub('^0[.]', '', thrs), NAME_LCM_ERR_POPUL_FILE_2, sep=''))
	# "lcm.err"

	foo.real = lcm.err[[j]]
	tmp.lambda = mean(foo.real)

	lambda.real = c(lambda.real, tmp.lambda)
      }
      rm(j)

      Lambda[[as.character(thrs)]]$real = lambda.real

      # >>> For the random data <<< #
      load(paste(NAME_LCM_RND_POPUL_DIR, NAME_LCM_RND_POPUL_FILE_0, 
	sub('^0[.]', '', thrs), NAME_LCM_RND_POPUL_FILE_2, sep=''))
      # [1] "lcm.rnd"

      lambda.rnd = sapply(lcm.rnd, mean)

      Lambda[[as.character(thrs)]]$rnd = lambda.rnd

      # >>> For Z-scores and P-values <<< #
      rnd.mean = mean(lambda.rnd)
      rnd.sd   = sd(lambda.rnd)

      z.scores = sapply(lambda.real, function(x.lambda) (x.lambda - rnd.mean) / rnd.sd)
      p.values = sapply(lambda.real, function(x.lambda) sum(lambda.rnd >= x.lambda) / length(lambda.rnd))

      results = cbind(Real_Err = lambda.real, Mean_Rnd = rnd.mean, SD_Rnd = rnd.sd, 
	Z_score = z.scores, P_value = p.values)

      Lambda[[as.character(thrs)]]$res = results
   }

   foo.lam = NULL
   for ( k in 1:length(Lambda) ) {
	tmp = cbind(as.numeric(names(Lambda)[k]), round(Lambda[[k]]$res, 3))
	foo.lam = rbind(foo.lam, tmp)
   }
   rm(k)
   colnames(foo.lam) = COLNAMES_ERR

   write.table(foo.lam, file = NAME_LAMBDA_ERR_RESULT_ZSC_PVAL_FILE[i], quote=F, row.names=F, sep='\t')
   save(Lambda, file = NAME_LAMBDA_ERR_RESULT_ALL_FILE[i])
}
rm(i)

cat('... End of PART XII ...\n\n')
setwd(CURRENT_DIR)


################################################################################

end.time = date()

cat('\n')
cat('> Program Start Time :', start.time, '\n')
cat('> Program End   Time :', end.time, '\n\n')
cat('Done!!!\n\n')


