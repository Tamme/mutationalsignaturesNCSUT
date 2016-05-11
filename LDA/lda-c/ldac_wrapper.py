import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import sys
sys.path.insert(1, '../')
sys.path.insert(1, '/home/tamme/Desktop/Masters/data')

#import file
from similarity import *
from scipy.spatial.distance import pdist, squareform


mat = scipy.io.loadmat('21breastWallHallsigs_4')
origMat = scipy.io.loadmat('21_WTSI_BRCA_whole_genome_substitutions.mat')
origMat2 = scipy.io.loadmat('21_WTSI_BRCA_whole_genome_substitutionsMutationOrder.mat')

#print mat
Wall = mat['Wall']
Hall = mat['Hall']
print Wall.shape
print Hall.shape
frob, avgStab, processStabAvg, H, exposureStd, W, centroidStd = evaluateStability(4, 100, Wall, Hall)

print processStabAvg
print H.shape
print W.shape
original = origMat['originalGenomes']
original2 = origMat2['originalGenomes']
print original.shape
difference =  W.dot(H) - original
difference2 = W.dot(H) - original2
		#print 'Error diff', difference.shape, np.sum(difference)
#frobenius = sqrt(sum(diag(t(diff) % * % diff)))
frob2 = np.sqrt(np.sum(difference.T.dot(difference).diagonal()))
frob3 = np.sqrt(np.sum(difference2.T.dot(difference2).diagonal()))
print frob2
print frob3
#scipy.io.savemat('WTEST.mat', mdict={'W': W, 'R':original, 'H':H})
#mutMatrix = mat['originalGenomes']
#mutMatrix = mutMatrix.T