#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import linalg as LA
from sklearn.preprocessing import normalize

ndim =36
oneday =8.64

def main():
#	C = np.arange(1,28,3).reshape(3,3) 
#	C = C.astype(np.float)
	C= np.random.normal(0.0,1.0,(ndim,ndim))
	dummy,C=LA.qr(C)
	print(C[:,5])
	C_sklearn=normalize(C,norm='l2',axis=0)
	C_coded=normalize_column(C)
	print("sklearn")
	print(C_sklearn[:,5])
	print("coded")
	print(C_coded[:,5])


def normalize_column(m):
	for i in range(ndim):
		print(i)
		c = m[:,i]
		print(c)
		lc = np.dot(c,c)**0.5
		m[:,i]/=lc
	return m
if __name__ == "__main__":
	main()
