import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import fetch_mldata
mnist = fetch_mldata('MNIST original')
X = mnist['data'] / 255.0
Xs = X

from rfn import *
W, P, Wout = train_rfn(Xs, 64, 500, 0.1, 0.1, 1e-1, 0.0, gpu_id="cpu")
print np.amin(W)

# plot weights
fig, ax = plt.subplots(5, 5, figsize=(8, 8))
for i, a in enumerate(ax.flat):
	a.pcolorfast(W[i].reshape(28, 28), cmap=plt.cm.Greys_r)
	a.set_ylim(28, 0)
	a.grid("off")
	a.set_axis_off()
fig.subplots_adjust(0, 0, 1, 1, 0, 0)
fig.show()
raw_input("Press Enter to continue...")

# calculate hidden units and reconstructions
H = np.maximum(0, np.dot(Wout, Xs.T))
R = np.dot(H.T, W)

# plot reconstructions
np.random.shuffle(R)  # shuffle samples, otherwhise we only plot 0s
fig, ax = plt.subplots(5, 5, figsize=(8, 8))
for i, a in enumerate(ax.flat):
	a.pcolorfast(R[i].reshape(28, 28), cmap=plt.cm.Greys_r)
	a.set_ylim(28, 0)
	a.grid("off")
	a.set_axis_off()
fig.subplots_adjust(0, 0, 1, 1, 0, 0)
fig.show()
raw_input("Press Enter to continue...")
np.save('W_mnist_out.npy', W)
np.save('H_mnist_out.npy', H)
np.save('R_mnist_out.npy', R)
