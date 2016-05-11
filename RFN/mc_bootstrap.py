import numpy as np

def monte_carlo_bootsrap(data_fs):
	#shape features x samples
	dshape = data_fs.shape
	mc_data = np.zeros(dshape)
	#print mc_data.shape
	for sample in range(dshape[0]):
		n = np.sum(data_fs[sample,])
		#print n
		prob = data_fs[sample,] / float(n)
		#print prob
		res = np.random.multinomial(n, prob, size=1)
		#print data_fs[sample,]
		#print res
		#print np.sum(res), np.sum(data_fs[sample,])
		#raw_input()
		mc_data[sample,] = res
	#print np.sum(mc_data)
	return mc_data

def removeWeak(data_fs):
	totalSum = np.sum(data_fs)
	removeWeak = 0.01 * totalSum
	previous = 0
	sums = np.sort(np.sum(data_fs, axis=1))
	idxes = np.argsort(np.sum(data_fs, axis=1))
	ctr = 0
	while True:
		previous += sums[ctr]
		if previous > removeWeak:
			break
		ctr += 1
	removeIdxes = idxes[:ctr]

	#print sums
	#print idxes
	#print removeIdxes
	data_fs = np.delete(data_fs, removeIdxes, axis=0)
	#print np.sum(data_fs)
	#print data_fs.shape
	return removeIdxes, data_fs

def addWeak(removeIdxes, W, Wstd, Wall, genomeErrors, genomeReconstructed):
	totalMutations = len(removeIdxes) + Wall.shape[0]
	W_final = np.zeros((totalMutations, W.shape[1]))
	Wstd_final = np.zeros((totalMutations, Wstd.shape[1]))
	Wall_final = np.zeros((totalMutations, Wall.shape[1]))
	genomeErrors_final = np.zeros((totalMutations, genomeErrors.shape[1], genomeErrors.shape[2]))
	genomeReconstructed_final = np.zeros((totalMutations, genomeReconstructed.shape[1], genomeReconstructed.shape[2]))
	origIdxCtr = 0
	for add in range(totalMutations):
		if add not in removeIdxes:
			W_final[add, :] = W[origIdxCtr, :]
			Wstd_final[add, :] = Wstd[origIdxCtr, :]
			Wall_final[add, :] = Wall[origIdxCtr, :]
			genomeErrors_final = genomeErrors[origIdxCtr, :, :]
			genomeReconstructed_final = genomeReconstructed[origIdxCtr, :, :]
			origIdxCtr += 1
	#print W_final

	return W_final, Wstd_final, Wall_final, genomeErrors_final, genomeReconstructed_final

def addWeakOnlyW(removeIdxes, W):
	totalMutations = len(removeIdxes) + W.shape[0]
	W_final = np.zeros((totalMutations, W.shape[1]))
	#Wstd_final = np.zeros((totalMutations, Wstd.shape[1]))
	#Wall_final = np.zeros((totalMutations, Wall.shape[1]))
	#genomeErrors_final = np.zeros((totalMutations, genomeErrors.shape[1], genomeErrors.shape[2]))
	#genomeReconstructed_final = np.zeros((totalMutations, genomeReconstructed.shape[1], genomeReconstructed.shape[2]))
	origIdxCtr = 0
	for add in range(totalMutations):
		if add not in removeIdxes:
			W_final[add, :] = W[origIdxCtr, :]
			#Wstd_final[add, :] = Wstd[origIdxCtr, :]
			#Wall_final[add, :] = Wall[origIdxCtr, :]
			#genomeErrors_final = genomeErrors[origIdxCtr, :, :]
			#genomeReconstructed_final = genomeReconstructed[origIdxCtr, :, :]
			origIdxCtr += 1
	#print W_final
	return W_final#, Wstd_final, Wall_final, genomeErrors_final, genomeReconstructed_final


#import scipy.io
#mat = scipy.io.loadmat('Breast/Breast_genomes_mutational_catalog_96_subs.mat')
#addWeak(mat['originalGenomes'])
