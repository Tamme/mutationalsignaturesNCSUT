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
from similarity import *
from rfn_transpose_many_iterations import *
import matplotlib.pyplot as plt
import pandas as pd
#http://pandas.pydata.org/pandas-docs/stable/visualization.html

if len(sys.argv) != 2:
	print('need nr of copies(python rfn_identify_nr_similarity.py nrOfCopies)')
	sys.exit()

#nrOfSignatures = int(sys.argv[1])
nrOfCopies = int(sys.argv[1])

frobs = []
stabs = []
sigRange = range(1, 20)
mat = scipy.io.loadmat('Breast/Breast_genomes_mutational_catalog_96_subs.mat')
mutMatrix_fs = mat['originalGenomes']
for nrOfSignatures in sigRange:	
	print nrOfSignatures
	frob, avgStab, _, _, _, _, _, _ = rfn_iterations_stability(mutMatrix_fs, nrOfSignatures, nrOfCopies)
	frobs.append(frob)
	stabs.append(avgStab)

#d = {'Frobenius' : frobs,
#     'Stability' : stabs
#     }
print frobs
print stabs
#frobs = [22748.3730882960226, 2268.5824501312832, 1754.2871049219652, 1470.5087237647117, 4256.6858378982597, 2109.9033911308502, 2871.0962150894102, 1451.323772353036]
#stabs = [0.99975652158012251, 0.9962228184689691, 0.97655222780720685, 0.99152053772610793, 0.72136133542582892, 0.75090715965392196, 0.61747043029221516, 0.79040179181477355]
#[24689.031394107522, 2746.9823148629657, 2229.1301771577387, 1803.626754449826, 1639.9129479531953, 3761.6297795244245, 3003.7740141786294, 2303.2489580446031, 2321.7273614689875, 1817.6064917843601, 2471.3823007599271, 3383.126110231759, 3808.3204811496703, 3613.6931796688564, 5060.4281512735752, 7570.6344825781562, 3129.080018721178, 5533.532014345511, 4024.6965951308889]
#[0.99996902763971929, 0.99987577094811142, 0.99469227436199492, 0.99397721380460102, 0.97253919670071343, 0.7020022121475189, 0.75267348277007085, 0.52271811372527821, 0.55906035498637097, 0.60495510850002454, 0.52184459676201567, 0.26516930070553063, 0.19520613716956198, 0.154744325456142, 0.046867955840426405, 0.11748658960798915, 0.16515625297231523, 0.1132471581310775, 0.064532361974746641]

#frobs = [2747.3153478633167, 2209.7375129084089]
#d = {'Stability': stabs, 'Frobenius': frobs}
#print d


#df = pd.DataFrame(d, index=list(sigRange))
#print
#df.plot(secondary_y=['Stability'], marker='o', color=['blue', 'black'])
#ax.set_xticks(list(sigRange))
#ax.set_ylabel('Frobenius')
#ax.right_ax.set_ylabel('Stability')
#ax.spines['left'].set_color('blue')
#ax.spines['right'].set_color('black')
#ax.yaxis.label.set_color('blue')
#plt.legend(loc='best')


f, axarr = plt.subplots(2, sharex=True)
axarr[0].set_ylim(0, 3000)
axarr[0].plot(list(range(1,len(frobs)+1)), frobs)
axarr[0].set_title('Frobenius error')
axarr[1].plot(list(range(1,len(frobs)+1)),stabs)
axarr[1].set_title('Stability')
axarr[0].xaxis.set_ticks(list(range(len(frobs)+2)))
axarr[1].xaxis.set_ticks(list(range(len(frobs)+2)))
plt.show()