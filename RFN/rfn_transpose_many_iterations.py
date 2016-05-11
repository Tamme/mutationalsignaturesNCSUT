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


def rfn_iterations_stability(mutMatrix_fs, nrOfSignatures, nrOfOutputs):
	nrOfIterations = 10000
	nrOfSamples = mutMatrix_fs.shape[1]
	saveOutput = False

	#nrOfSignatures = 10
	#nrOfOutputs = 8

	#        _librfn.train_cpu(X, W, P, X.shape[0], X.shape[1], n_hidden, n_iter,
	#                          batch_size, etaW, etaP, minP, h_threshold, dropout_rate, input_noise_rate,
	#                          l2_weightdecay, l1_weightdecay, momentum, _input_noise_types[input_noise_type],
	#                          _activation_types[activation], 1, applyNewtonUpdate, seed)
	mutMatrix_fs_orig = mutMatrix_fs.copy()
	removeIdxes, mutMatrix_fs = removeWeak(mutMatrix_fs)
	nrOfMutations = mutMatrix_fs.shape[0]

	Wall  = np.zeros((nrOfMutations, nrOfSignatures*nrOfOutputs))
	Hall  = np.zeros((nrOfSignatures*nrOfOutputs, nrOfSamples))
	genomeErrors = np.zeros((nrOfMutations, nrOfSamples, nrOfOutputs))
	reconstructed = np.zeros((nrOfMutations, nrOfSamples, nrOfOutputs))

	for outputNr in range(nrOfOutputs):
		mutMatrixIter = monte_carlo_bootsrap(mutMatrix_fs)
		results = ''
		minF = -1
		maxF = -1

		for x in range(15):
			W, P, Wout = train_rfn(mutMatrixIter, nrOfSignatures, nrOfIterations, 0.1, 0.1, 1e-1, 0.0, l1_weightdecay=0.00001, gpu_id="cpu")


			H = np.maximum(0, np.dot(Wout, mutMatrixIter.T))
			R = np.dot(H.T, W)

			difference =  np.abs(mutMatrixIter - R)
			#print 'Error diff', difference.shape, np.sum(difference)
			frob = np.sqrt(np.sum((difference.T.dot(difference)).diagonal()))
			#print 'Frobenius \t\t\t\t\t\t\t',  frob
			results += str(frob) + '  '
			if minF == -1 or minF > frob:
				minF = frob
				minR = R
				minH = H
				minW = W
			if maxF == -1 or maxF < frob:
				maxF = frob
				maxR = R
		#print results
		#print H.shape
		print outputNr, 'min', minF, 'max', maxF
		H = minH
		W = minW
		R = minR
		
		#print 'H', H.shape, np.amin(H)
		#change W and H back from transposing
		tmp = W
		W = H
		H = tmp

		#print W.shape
		print 'nr of sparse', np.sum(W == 0.0), np.sum(W == 0.0)/float(W.shape[0]*W.shape[1])
		

		reconstructed[:, :, outputNr] = R
		genomeErrors[:, :, outputNr] = difference
		Wall[:, (outputNr*nrOfSignatures):((outputNr*nrOfSignatures)+nrOfSignatures)] = W.T
		Hall[(outputNr*nrOfSignatures):((outputNr*nrOfSignatures)+nrOfSignatures), :] = H


	print 'Final'

	frob, avgStab, processStabAvg, H, Hstd, W, Wstd = evaluateStability(nrOfSignatures, nrOfOutputs, Wall, Hall)
	W, Wstd, Wall, genomeErrors, genomeReconstructed = addWeak(removeIdxes, W, Wstd, Wall, genomeErrors, reconstructed)
	R = np.dot(W, H)
	difference = mutMatrix_fs - R
	frob = np.sqrt(np.sum((difference.T.dot(difference)).diagonal()))
	scipy.io.savemat('to_matlab/All_breast_MC_stability_C_L1_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_' + str(nrOfOutputs) + '_out.mat', 
	mdict={'W': W, 'H':H, 'Wstd': Wstd, 'Hstd': Hstd, 'stability':processStabAvg, 'Wall':Wall, 'Hall':Hall})
	
	return frob, avgStab, processStabAvg, H, Hstd, W, Wstd, R


if __name__ == "__main__":
	mat = scipy.io.loadmat('Breast/Breast_genomes_mutational_catalog_96_subs.mat')
	mutMatrix_fs = mat['originalGenomes']
	#mutMatrix = mutMatrix.T
	#print mutMatrix.shape
	#raw_input()

	if len(sys.argv) != 3:
		print('need nr of signatures and ierations(python my_rfn....py nrOfsignatures nrOfOutputs)')
		sys.exit()

	nrOfSignatures = int(sys.argv[1])
	nrOfOutputs = int(sys.argv[2])

	frob, avgStab, processStabAvg, H, Hstd, W, Wstd, R = rfn_iterations_stability(mutMatrix_fs, nrOfSignatures, nrOfOutputs)

	print processStabAvg
	print frob
	#raw_input()
	#if saveOutput == True:
	#scipy.io.savemat('to_matlab/All_breast_MC_C_L1_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_' + str(nrOfOutputs) + '_out.mat', 
	#%mdict={'Wall': Wall, 'Hall':Hall, 'reconstruct':reconstructed, 'genomeErrors':genomeErrors})
	#np.save('models/W_breast_C_L1_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_' + str(nrOfOutputs) + '_out.npy', W)
	#scipy.io.savemat('models/breast_C_L1_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_' + str(nrOfOutputs) + '_out.mat', mdict={'W': W, 'H':H, 'R':R})
	raw_input("Press Enter to continue...")
