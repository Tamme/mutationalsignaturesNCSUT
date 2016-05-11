#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from cpu_final_basic_python_implementation import *


mat = scipy.io.loadmat('21_WTSI_BRCA_whole_genome_substitutions.mat')
#print 
mutMatrix = mat['originalGenomes'] 
mutMatrix = mutMatrix.T

mat = scipy.io.loadmat('toy_data_3.mat')
mutMatrix = mat['X']

#normalize
#print mutMatrix.shape
#mutMatrix = mutM

print mutMatrix.shape
mutMatrix = mutMatrix / np.sum(mutMatrix, axis=0)

print np.sum(mutMatrix, axis=0)

#centering didnt help
#mean = mutMatrix.mean(axis=1) 
#mutMatrix = mutMatrix - mean[:, np.newaxis]


#from sklearn.datasets import fetch_mldata
#mnist = fetch_mldata('MNIST original')
#X = mnist['data'] / 255.0
#print 'X', X.shape
#print X.shape
#print X[1:30, 1:20]
#print sum(X)

#Plot data
#nrOfSamples = 21
#f, axarr = plt.subplots(nrofSamples, sharex=True)
#for i in range(nrofSignatures):
#	axarr[i].bar(range(96), mutMatrix[i, :])
#fig = plt.gcf()
#fig.show()
#raw_input("Press Enter to continue...")

#W, P, Wout = train_rfn(mutMatrix, nrofSignatures, 500, 0.1, 0.1, 1e-1, 0.0, gpu_id="cpu")
nrOfIterations = 10000
nrOfSignatures = 7
nrOfSamples = 500

#lol = mutMatrix[:-6,:]
#mutMatrix = np.delete(mutMatrix, 15, 0)
#print mutMatrix.shape
#print np.sum(mutMatrix, axis=0)
#print np.sum(mutMatrix, axis=1)



W, H, P = train_rfn_cpu(mutMatrix.T, nrOfSignatures, nrOfIterations, 0.1, 0.1, 0.0)
#W, P, Wout = train_rfn(mutMatrix, nrofSignatures, nrOfIterations, 0.1, 0.1, 1e-1, 0.0, gpu_id="cpu")


#H = np.maximum(0, np.dot(Wout, mutMatrix.T))
R = np.dot(H.T, W)

print 'W', W.shape, np.amin(W)
print 'H', H.shape, np.amin(H)
print 'R', R.shape, np.amin(R)
print 'data', mutMatrix.shape

print '-------------------'

#print H 

difference =  np.abs(mutMatrix.T - R)
print 'Error diff', np.sum(difference)
print 'Frobenius norm',  np.sqrt(np.sum((difference.T.dot(difference)).diagonal()))


#print H
print H.shape
f, axarr = plt.subplots(nrOfSignatures, sharex=True)
for i in range(nrOfSignatures):
	axarr[i].bar(range(96), H[i, :])
fig = plt.gcf()
#fig.show()

print W.shape
f, axarr = plt.subplots(nrOfSignatures, sharex=True)
for i in range(nrOfSignatures):
	axarr[i].bar(range(nrOfSamples), W[i, :])
fig = plt.gcf()
#fig.show()

#plt.bar(range(96), W[1, :])
#print 'max', max(W[1,:])
#print 'min', min(W[1,:])

#fig
#print P.shape
#print Wout.shape
#print W[:,70:80]
#return
# plot weights

np.save('W_breast_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_out.npy', W)
np.save('H_breast_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_out.npy', H)
np.save('R_breast_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_out.npy', R)
scipy.io.savemat('W_breast_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_out.mat', mdict={'W': W})
scipy.io.savemat('H_breast_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_out.mat', mdict={'H': H})
scipy.io.savemat('R_breast_' + str(nrOfSignatures) + '_' + str(nrOfIterations) + '_out.mat', mdict={'R': R})
raw_input("Press Enter to continue...")
