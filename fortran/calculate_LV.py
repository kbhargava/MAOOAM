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
#	traj, tlm, tim = read_file("evol_field_tlm.dat")
	dt=1
	tmax=300000
#	dt = tim[1]-tim[0]
	print(dt)
	step = np.int(dt)
#	tmax = len(tim)
	t_0 = np.int(tmax/3)
	t_N = 2*t_0
	t_w = tmax
	t_c = tmax-t_0
	rescale = 1#np.int(1*dt)
	
#	lyap_sum =np.zeros(ndim)
#	cnt = 0
#	M = np.eye(ndim)
#	S = np.eye(ndim)
#	LE = np.empty(ndim)
#	Qhist = np.empty((tmax, ndim, ndim))
#	Rhist = np.empty((tmax, ndim, ndim))
#	Chist = np.empty((tmax, ndim, ndim))
	CLV = np.empty((t_N-t_0, ndim, ndim))
	Mhist = np.empty((t_N-t_0, ndim, ndim))
#
## Forward Integration	
#	for i in range(0,tmax):
#		M = np.dot((tlm[i,:,:]),M)
#		if (np.mod(i+1,rescale)==0):
#			Q,R = LA.qr(M)
#			Qhist[i,:,:] = Q
#			Rhist[i,:,:] = R
#			M = Q
#			rdiag = abs(np.diag(R))
#			lyap_sum =lyap_sum+np.log(rdiag)/(rescale/oneday)
#			cnt = cnt+1
#	print(cnt)
#
## Computing Lyapunov Exponents	 
#	LE = lyap_sum/cnt
##	print (LE)
##	np.savetxt('LE_120000.txt',LE)
##	plt.plot(LE)
##	plt.savefig("LE120000.png",dpi = 300)
#
## Initializing C matrix
#	C = np.random.normal(0.0,1.0,(ndim,ndim))
#	C = np.triu(C)
#
## Backward Integration
#	for i in range(tmax-1,t_0, -1):
##		Chist[i,:,:]=normalize(C,norm='l2',axis=0)
#		Chist[i,:,:]=normalize_column(C)
#		R_inv = LA.inv(Rhist[i-1,:,:])
#		C = np.dot(R_inv,C)
#	Chist[t_0,:,:]=normalize_column(C)
#	
## Calculating Lyapunov Vectors
#	for idx,i in enumerate(range(t_0,t_N,1)):
#		print(idx,i)
#		Q = Qhist[i,:,:]
#		C = Chist[i,:,:]
#		CLV[idx,:,:] = np.dot(Q,C)
#		Mhist[idx,:,:] = tlm[i,:,:]
#
#
## Saving data
#	np.save("CLV.npy",CLV)
#	np.save("Mhist.npy",Mhist)

	CLV=np.load("CLV.npy")
	Mhist=np.load("Mhist.npy")
	print("Loaded")
# Calculating growth rates of the error along different CLVs
	CLV_list=[0,1,4,9,14,19,24,29] #CLV number to compute growth rate for
	del_t = 100 # time to compute growth rates for
	norm_vi = np.empty((t_N-t_0-del_t,del_t+1,len(CLV_list)))
	for idx,CLV_num in enumerate(CLV_list):
		for t_index in range(t_N-t_0-del_t):
			print(t_index)
			vi = CLV[t_index,:,CLV_num]
			norm_vi[t_index,0,idx] = LA.norm(vi)
			for dtx in range(del_t):
				vi = np.dot(Mhist[t_index+dtx,:,:],vi[:])
				norm_vi[t_index,dtx+1,idx] = LA.norm(vi)
	print(norm_vi[t_N-t_0-del_t-1,:,0])
	norm_log_mean = np.mean(np.log(norm_vi),axis = 0)
	norm_mean_log = np.log(np.mean(norm_vi,axis=0))
	growth_log_mean = np.copy(norm_log_mean)
	growth_mean_log = np.copy(norm_mean_log)
	for dtx in reversed(range(del_t)):
		growth_log_mean[dtx+1,:] -= growth_log_mean[dtx,:]
		growth_mean_log[dtx+1,:] -= growth_mean_log[dtx,:]
	

# Plotting growth rates of the error along different CLVs
	x = (np.arange(float(del_t)+0.5))*(dt/oneday)
	for idx,CLV_num in enumerate(CLV_list):
		plt.plot(x,growth_log_mean[1:,idx], label="CLV num%d" %(CLV_num+1))
	plt.xlabel("days")
	plt.legend()
	plt.savefig("fig8_growth_log_mean.png",dpi=300)
	plt.close()


	for idx,CLV_num in enumerate(CLV_list):
		plt.plot(x,growth_mean_log[1:,idx], label="CLV num%d" %(CLV_num+1))
	plt.xlabel("days")
	plt.legend()
	plt.savefig("fig8_growth_mean_log.png",dpi=300)
	plt.close()


def normalize_column(m):
	for i in range(ndim):
		c = m[:,i]
		lc = np.dot(c,c)**0.5
		m[:,i]/=lc
		return m

def read_file(file):
	#return np.ndarray[time, ndim]
	with open(file, "r") as f:
		ar = f.read().split()
	ar2 = []
	n = len(ar)
	nrec = ndim ** 2 + ndim + 1
	tmax2=n/nrec
	print (tmax2)
	na = np.empty((tmax2, ndim))
	ntl = np.empty((tmax2, ndim ** 2))
	tim = np.empty((tmax2))
	
	for i in range(tmax2):
		tim[i]   = ar[i * nrec]
		na[i, :] = ar[i * nrec + 1:i * nrec + ndim + 1]
		ntl[i, :] = ar[i * nrec + ndim + 1:i * nrec + ndim ** 2 + ndim + 1]
	ntl = ntl.reshape((tmax2, ndim, ndim))
	return na, ntl, tim

if __name__ == "__main__":
	main()
