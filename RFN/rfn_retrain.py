#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from rfn import *
from mc_bootstrap import *
from similarity import *
from scipy.spatial.distance import pdist, squareform
import copy

#mat = scipy.io.loadmat('21_WTSI_BRCA_whole_genome_substitutions.mat')


def rfn_iterations_stability(mutMatrix_fs, nrOfSignatures, nrOfOutputs, saveOutput = False, l1 = 0.000001):
	nrOfIterations = 100 #important, small one didnt make it converge
	nrOfSamples = mutMatrix_fs.shape[1]
	 #= False

	#nrOfSignatures = 10
	#nrOfOutputs = 8

	#        _librfn.train_cpu(X, W, P, X.shape[0], X.shape[1], n_hidden, n_iter,
	#                          batch_size, etaW, etaP, minP, h_threshold, dropout_rate, input_noise_rate,
	#                          l2_weightdecay, l1_weightdecay, momentum, _input_noise_types[input_noise_type],
	#                          _activation_types[activation], 1, applyNewtonUpdate, seed)
	mutMatrix_fs_orig = mutMatrix_fs.copy()
	removeIdxes, mutMatrix_fs = removeWeak(mutMatrix_fs)
	nrOfMutations = mutMatrix_fs.shape[0]

	bestF = -1
	for getBest in range(1, 10):
		for outputNr in range(nrOfOutputs):
			mutMatrixIter = monte_carlo_bootsrap(mutMatrix_fs)
			results = ''
			minF = -1
			maxF = -1
			if outputNr == 0:
				nonConvIter = 200 #need to search a lot initially
			else:
				nonConvIter = 1


			for x in range(nonConvIter):
				if outputNr != 0:
					W, P, Wout = train_rfn(mutMatrixIter, nrOfSignatures, nrOfIterations, 0.01, 0.01, 1e-6, 0.0, startW=Wold, l1_weightdecay=l1, gpu_id="cpu")			
				else:
					W, P, Wout = train_rfn(mutMatrixIter, nrOfSignatures, nrOfIterations, 0.1, 0.1, 1e-2, 0.0, l1_weightdecay=l1, gpu_id="cpu")
				H = np.maximum(0, np.dot(Wout, mutMatrixIter.T))
			
				R = np.dot(H.T, W)

				difference =  np.abs(mutMatrixIter - R)
				frob = np.sqrt(np.sum((difference.T.dot(difference)).diagonal()))
				results += str(frob) + '  '
				if minF == -1 or minF > frob:
					minF = frob
					minR = R
					minH = H
					minW = W
				if maxF == -1 or maxF < frob:
					maxF = frob
					maxR = R
			H = minH
			W = minW
			Wold = minW
			R = minR
			
		if bestF == -1 or bestF > minF:
			bestH = H
			bestW = W
			bestF = minF

	print 'Final'
	tmp = bestW
	W = bestH.T
	H = tmp
	#frob, avgStab, processStabAvg, H, Hstd, W, Wstd = evaluateStability(nrOfSignatures, nrOfOutputs, Wall, Hall)
	W = addWeakOnlyW(removeIdxes, W)
	R = np.dot(W, H)
	difference = mutMatrix_fs_orig - R
	frob = np.sqrt(np.sum((difference.T.dot(difference)).diagonal()))
	print frob
	if saveOutput == True:
		scipy.io.savemat('to_matlab/All_breast_MC_retrain_C_L1_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_' + str(nrOfOutputs) + '_out.mat', mdict={'W': W, 'H':H, 'R':R})
	
	return frob, W, H, R


if __name__ == "__main__":
	mat = scipy.io.loadmat('Breast/Breast_genomes_mutational_catalog_96_subs.mat')
	#mat = scipy.io.loadmat('21_WTSI_BRCA_whole_genome_substitutions.mat')
	#mat = scipy.io.loadmat('506large96SubsData.mat')
	
	mutMatrix_fs = mat['originalGenomes']
	#mutMatrix = mutMatrix.T
	#print mutMatrix.shape
	#raw_input()

	if len(sys.argv) != 3:
		print('need nr of signatures and ierations(python my_rfn....py nrOfsignatures nrOfOutputs)')
		sys.exit()

	nrOfSignatures = int(sys.argv[1])
	nrOfOutputs = int(sys.argv[2])

	frob, W, H, R = rfn_iterations_stability(mutMatrix_fs, nrOfSignatures, nrOfOutputs)

	#raw_input()
	#if saveOutput == True:
	#scipy.io.savemat('to_matlab/All_breast_MC_C_L1_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_' + str(nrOfOutputs) + '_out.mat', 
	#%mdict={'Wall': Wall, 'Hall':Hall, 'reconstruct':reconstructed, 'genomeErrors':genomeErrors})
	#np.save('models/W_breast_C_L1_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_' + str(nrOfOutputs) + '_out.npy', W)
	#scipy.io.savemat('models/breast_C_L1_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_' + str(nrOfOutputs) + '_out.mat', mdict={'W': W, 'H':H, 'R':R})
	raw_input("Press Enter to continue...")
