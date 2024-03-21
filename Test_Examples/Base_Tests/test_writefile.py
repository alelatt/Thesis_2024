import numpy as np
import time as tm
import os

t1 = tm.time()
for i in range(0,100):
	times = np.arange(0, int(1e4), 0.1)
	semi = []
	ecc = []
	inc = []

	for i in range(len(times)):
		semi.append([1.2 + i/2, 2.2 + i/2, 3.2 + i/2])
		ecc.append([0 + i/3, 0 + i/2, 0 + i])
		inc.append([0, 90, 180])

	out = [times, semi, ecc, inc]

	if os.path.isfile("./out.npz"):
		outs = np.load("./out.npz")
		np.savez_compressed("./out.npz", t = np.concatenate((outs['t'], np.array(out[0]))), a = np.concatenate((outs['a'], np.array(out[1]))), e = np.concatenate((outs['e'], np.array(out[2]))), i = np.concatenate((outs['i'], np.array(out[3]))))
	else:
		np.savez_compressed("./out.npz", t = np.array(out[0]), a = np.array(out[1]), e = np.array(out[2]), i = np.array(out[3]))

print(tm.time()-t1)

outs = np.load("./out.npz")
print(len(outs['t']))