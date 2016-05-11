from my_rfn_retrain import *
import numpy as np

if __name__ == "__main__":
	#mat = scipy.io.loadmat('Breast/Breast_genomes_mutational_catalog_96_subs.mat')
	#mat = scipy.io.loadmat('21_WTSI_BRCA_whole_genome_substitutions.mat')
	mat = scipy.io.loadmat('506large96SubsData.mat')
	
	mutMatrix_fs = mat['originalGenomes']
	#mutMatrix = mutMatrix.T
	#print mutMatrix.shape
	#raw_input()

	if len(sys.argv) != 2:
		print('need nr of iterations(python my_rfn....py nrOfOutputs)')
		sys.exit()

	nrOfOutputs = int(sys.argv[1])

	for signature in [7]:
		print signature, '----------'
		sparsities = np.arange(0, 5.1, 0.05)#(0.1, )list(range(1, 101)) / float(10)
		#print sparsities
		#print 9 / float(10)
		#raw_input()
		for sparsity in sparsities:
			print sparsity
			frob, W, H, R = rfn_iterations_stability(mutMatrix_fs, signature, nrOfOutputs, False, sparsity)

		#print frob
		
			scipy.io.savemat('models/506_breast_C_L1_S_' + str(sparsity) + '_' + str(signature) + '_' + str(nrOfOutputs) + '_out.mat', mdict={'W': W, 'H':H, 'R':R})
	raw_input("Press Enter to continue...")