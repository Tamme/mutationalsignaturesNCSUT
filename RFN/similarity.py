import numpy as np
from scipy.spatial.distance import squareform, pdist
from sklearn.metrics import silhouette_samples
from scipy.cluster.vq import kmeans2
from scipy.spatial.distance import pdist, squareform
from sklearn import datasets
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.io import loadmat
import copy
import sys
import scipy.io


def silhouette(X, cIDX):
    """
    Computes the silhouette score for each instance of a clustered dataset,
    which is defined as:
        s(i) = (b(i)-a(i)) / max{a(i),b(i)}
    with:
        -1 <= s(i) <= 1

    Args:
        X    : A M-by-N array of M observations in N dimensions
        cIDX : array of len M containing cluster indices (starting from zero)

    Returns:
        s    : silhouette value of each observation
    """

    N = X.shape[0]              # number of instances

    K = len(np.unique(cIDX))    # number of clusters

    # compute pairwise distance matrix
    D = squareform(pdist(X))
    print D
    print N
    print K
    #print cIDX[3]

    # indices belonging to each cluster
    kIndices = [np.flatnonzero(cIDX==k) for k in range(K)]

    # compute a,b,s for each instance
    a = np.zeros(N)
    b = np.zeros(N)
    print "len", len(kIndices)
    print kIndices
    for i in range(N):
        # instances in same cluster other than instance itself
#	print int(cIDX[i][0])
	#for ind in kIndices[int(cIDX[i][0])]:
#		if ind!=i:#
	#	print ind
#	print [D[i][ind] for ind in kIndices[cIDX[i]] if ind!=i]
        a[i] = np.mean( [D[i][ind] for ind in kIndices[int(cIDX[i][0])] if ind!=i] )
        # instances in other clusters, one cluster at a timed
        b[i] = np.min( [np.mean(D[i][ind]) 
		        for k,ind in enumerate(kIndices) if cIDX[i]!=k] )
    	s = (b-a)/np.maximum(a,b)
	return s
    print "0000000"

	




#nrOfSignatures = 10
#nrOfCopies = 8

def evaluateStability(nrOfSignatures, nrOfCopies, dataW = [], dataH = []):

	BIG_NUMBER = 100
	CONVERG_ITER = 10
	CONVERG_CUTOFF = 0.005
	TOTAL_INIT_CONDITIONS = 5
	totalProcesses = nrOfSignatures
	totalReplicates = 100 # 100
	processesDist = 'cosine'
	if len(dataW) == 0:
		mats = loadmat("to_matlab/All_breast_MC_C_L1_" + str(nrOfSignatures) + "_10000_" + str(nrOfCopies) + "_out.mat")
		Wall = mats["Wall"]
		Hall = mats["Hall"]
	else:
		Wall = dataW
		Hall = dataH


	for j in range (nrOfSignatures*nrOfCopies):
		#print j
		total = np.sum(Wall[:, j])
		add = 0.00000
		if total == 0.0:
			Wall[:, j] += add
			#add += 0.00001
			total = np.sum(Wall[:, j])
			#print j, total
		Wall[:,j] = Wall[:,j] / total
		Hall[j,:] = Hall[j,:] * total
	#print Wall.shape
	#print Hall.shape
	#print Wall
	#print Hall
	#raw_input()

	#print Hall
	#raw_input()

	#Wall = np.random.rand(96,nrOfSignatures*nrOfCopies)
	#Hall = np.random.rand(nrOfSignatures*nrOfCopies, 21)
	minClusterDist = BIG_NUMBER
	totalIter = Wall.shape[1] / nrOfSignatures
	idx = np.zeros((Hall.shape[0], 1))
	clusterCompactness = np.zeros((nrOfSignatures, totalIter))
	iStartigDataSet = list(range(0, Wall.shape[1], nrOfSignatures))
	#print iStartigDataSet	
	iStartingDataSet = np.random.permutation(iStartigDataSet)
	#print iStartigDataSet
	#print iStartingDataSet
	#raw_input()

	for iInitData in range(0, np.minimum(TOTAL_INIT_CONDITIONS, totalIter)):
		#assignin('base', 'iInitData', iInitData);
		iStartingData = iStartingDataSet[iInitData]
		#disp(iStartingData);
		iEnd = iStartingData + totalProcesses
		#print iStartingData, iEnd


		centroids = Wall[:, iStartingData:iEnd].copy()
		#print 'ce', centroids.shape
		#print centroids
		#raw_input()
		centroidsTest = np.ones((centroids.shape[0], centroids.shape[1]))
		countIRep = 0
		for iRep in range(totalReplicates):
			tmp = (np.concatenate((centroids, Wall), axis = 1)).T
			#print Wall
			#print tmp
			allDist = squareform(pdist(tmp, processesDist) );
			centroidDist = allDist[0:centroids.shape[1], (centroids.shape[1]):allDist.shape[1]].T
			#print centroidDist
			#print centroidDist.shape
			#raw_input()
			#print centroidDist.shape
			jRange = np.random.permutation(totalProcesses)
			for jIndex in range(totalProcesses):
				j = jRange[jIndex]
				for i in range(0, Wall.shape[1], totalProcesses):
					iRange = list(range(i, i+totalProcesses))
					#print 'ir'				
					#print iRange
					#print j, iRange
					#print centroidDist.shape
					#print centroidDist[iRange, j]

					Ind = np.argmin(centroidDist[iRange, j])
					#print 'Ind', Ind
					centroidDist[iRange[Ind],:] = BIG_NUMBER
					idx[iRange[Ind]] = j
			#print idx
			#raw_input()
			
			#print '////////////////////////////'
			maxDistToNewCentroids =  0
			#print centroids.shape[1]
			for i in range(centroids.shape[1]):
				same = np.where(idx == i)[0]
	#			print np.mean(Wall[:, same], axis = 1)
				#print idx == i
				#print Wall
				#print 'wall'
				#print Wall
				#print same
				#print 'same'
				centroids[:, i] = np.mean(Wall[:, same], axis = 1)
				#print centroids[:, i]
				tmp = (np.concatenate(([centroids[:, i, np.newaxis]], [centroidsTest[:, i, np.newaxis]]), axis = 0))
				
				#print tmp
				#raw_input()
				tmp = tmp.reshape((tmp.shape[0],tmp.shape[1])) #TODOOO
				#print tmp
				#print tmp.shape
				#raw_input()
				#r = pdist(tmp, 'cosine')
				#print r
				maxDistToNewCentroids = max(maxDistToNewCentroids, pdist(tmp, 'cosine'))
			
	#		end
			#print maxDistToNewCentroids
	#		print centroids
			#raw_input()
		
			if maxDistToNewCentroids < CONVERG_CUTOFF:
				countIRep += 1
			else:
				countIRep = 0
				centroidsTest = centroids.copy()
			#print countIRep
			#print centroids.shape
			#print centroids
			#raw_input()

			if countIRep == CONVERG_ITER:
				break

		#end
		for i in range(centroids.shape[1]):
			same = np.where(idx == i)[0]
	#		centroids
			#print centroids[:, i, np.newaxis].shape, Wall[:, same]
			tmp = np.concatenate((centroids[:, i, np.newaxis], Wall[:, same]), axis=1).T
			#print tmp.shape
			#print tmp
			clusterDist = squareform(pdist(tmp, 'cosine'))
			#print clusterDist
			#raw_input()
	 		clusterCompactness[i,:] = clusterDist[0, 1:clusterDist.shape[1]]
		#print clusterCompactness
		#print centroids

		#print idx
		#print Wall
		
		#raw_input()
		#end
		#print 'mean', np.mean(clusterCompactness)
		if minClusterDist > np.mean(clusterCompactness):
			centroidsFinal = centroids.copy()
			idxFinal = idx.copy()
			clusterCompactnessFinal = clusterCompactness.copy()
		#print 'ctrr'
	#	print centroids

	#end
	centroids = centroidsFinal.T
	#print centroids
	#print idxFinal

	idx = idxFinal.copy()
	clusterCompactness = clusterCompactnessFinal

	centDist = np.mean(clusterCompactness, axis=1)
	#print centDist
	#print 'ww'

	centDistInd = np.argsort(centDist)
	#print centDistInd
	#print idx
	#raw_input()
	clusterCompactness = clusterCompactness[centDistInd, :]
	#print clusterCompactness
	centroids = centroids[centDistInd, :]
	idxNew = idx.copy()
	#print 'idxNew', idxNew
	#print 'eeee'
	#print centDistInd
	#raw_input()

	for i in range(totalProcesses):##
		same = np.where(idx == centDistInd[i])[0]#
		idxNew[same] = i

	idx = idxNew
	#print idx
	#raw_input()
	#print idx
	#raw_input()

	#print idx
	if totalProcesses > 1: #TODOO ELSE
		idx = idx.flatten()
		#D = squareform(pdist(Wall.T, metric="cosine"))
		processStab = silhouette_samples(Wall.T, idx, metric='cosine')
		#processStab = np.round(processStab, 4)
		#print processStab
		processStabAvg = np.zeros((totalProcesses))
		#print processStabAvg.shape
		for i in range(totalProcesses):
			same = np.where(idx == i)[0]
			#print processStab[same]
			processStabAvg[i] = np.mean(processStab[same])
	#print 'Stability', processStabAvg
	else:
		same = np.where(idx == i)[0]
	#		centroids
			#print centroids[:, i, np.newaxis].shape, Wall[:, same]
		tmp = np.concatenate((centroids.T, Wall), axis=1).T
			#print tmp.shape
			#print tmp
		allDist = squareform(pdist(tmp, 'cosine'))
		processStab = 1 - allDist[0:centroids.T.shape[1], centroids.T.shape[1]:allDist.shape[1]].T
		processStabAvg = np.mean(processStab)
		#print processStab
		#raw_input()


	centroidStd = np.zeros((centroids.shape[0], centroids.shape[1],))
	for i in range(totalProcesses):
		same = np.where(idx == i)[0]
		centroidStd[i,:] = np.std(Wall[:, same], axis = 1, ddof=1)

	centroids = centroids.T
	centroidStd = centroidStd.T

	idxS = np.zeros((idx.shape[0], 1))
	for i in range(0, Wall.shape[1], totalProcesses):
		iEnd = i + totalProcesses
		idxG = idx[i:iEnd]
		for j in range(totalProcesses):
			idxS[i+j, :] = np.where(idxG == j)[0]
	#print iEnd

	#print 'idsxs', idxS
	#

	exposure = np.zeros((int(max(idxS.flatten()))+1, Hall.shape[1]))
	exposureStd = np.zeros((int(max(idxS.flatten()))+1, Hall.shape[1]))
	#print exposure.shape

	for i in range(int(max(idxS.flatten()))+1):
		same = np.where(idx == i)[0]
		#print same
		#print Hall[same, :]
		#print np.mean(Hall[same, :], axis = 0)
		#raw_input()
		exposure[i, :] = np.mean(Hall[same, :], axis = 0)

		exposureStd[i, :] = np.std(Hall[same, :], axis = 0, ddof=1)
	
	#print Hall


	#print exposure.shape
	#print centroids.shape

	#compare to the original
	#mat = loadmat('21_WTSI_BRCA_whole_genome_substitutions.mat')
	#mutMatrix = mat['originalGenomes']
	#mutMatrix = mutMatrix.T

	#R = np.dot(centroids, exposure)

	#difference =  np.abs(mutMatrix - R)
			#print 'Error diff', difference.shape, np.sum(difference)
	#frob = np.sqrt(np.sum((difference.T.dot(difference)).diagonal()))
	#print centroidStd
	#raw_input()
	frob = 0
	avgStab = np.mean(processStabAvg)
	#print frob, avgStab
	return frob, avgStab, processStabAvg, exposure, exposureStd, centroids, centroidStd

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print('need nr of signatures and ierations(python similarity.py nrOfsignatures nrOfOutputs)')
		sys.exit()

	nrOfSignatures = int(sys.argv[1])
	nrOfCopies = int(sys.argv[2])
	frob, avgStab, processStabAvg, H, W = evaluateStability(nrOfSignatures, nrOfCopies)
	print H.shape, W.shape
	R = np.dot(W, H)

	scipy.io.savemat('models/sim_res_breast_C_L1_' + str(nrOfSignatures) + '_' + str(nrOfCopies) + '_out.mat', mdict={'W': W, 'R':R, 'H':H})
	print processStabAvg
