from rfn_retrain import *

if __name__ == "__main__":
	#mat = scipy.io.loadmat('data/119_genomes_96_subs_data.mat')
	#mat = scipy.io.loadmat('506_genomes_96_subs_data.mat')
	#mat = scipy.io.loadmat('data/synthetic_5sigs_96mut_500samples_2000tot_20samp.mat')
	mat = scipy.io.loadmat('data/synthetic_5sigs_96mut_2000tot_20samp.mat')
	
	mutMatrix_fs = mat['originalGenomes']
	
	if len(sys.argv) != 2:
		print('need nr of signatures and ierations(python my_rfn....py nrOfOutputs)')
		sys.exit()

	nrOfOutputs = int(sys.argv[1])

	for signature in range(3,8):
		print signature, '----------'
		frob, W, H, R = rfn_iterations_stability(mutMatrix_fs, signature, nrOfOutputs)
		print frob
		scipy.io.savemat('models/2synthetic_2000_20_C_L1_0_' + str(signature) + '_' + str(nrOfOutputs) + '_out.mat', mdict={'W': W, 'H':H, 'R':R})
	raw_input("Press Enter to continue...")