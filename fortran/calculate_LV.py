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
#--------------------------------------------------------------------------------
#			Setting up variables
#--------------------------------------------------------------------------------
	traj, tlm, tim = read_file("evol_field_tlm_takuma_ic_2.dat")
#	traj, tlm, tim = read_file("evol_field_tlm_120000.dat")
	dt = tim[1]-tim[0]
	step = np.int(np.round(dt))
	tmax = len(tim)
	rescale = 10 
	t_res=np.int(np.round(tmax/rescale))
	t_0 = np.int(t_res/3)+1
	t_N = 2*(t_0-1)
	t_w = t_res
	t_c = t_res-t_0+1

	print("tmax,dt"+str(tmax)+","+str(rescale))
	lyap_sum =np.zeros(ndim)
	cnt = 0
	M = np.eye(ndim)
	S = np.eye(ndim)
	LE = np.empty(ndim)
	Qhist = np.empty((t_res, ndim, ndim))
	Rhist = np.empty((t_res, ndim, ndim))
	Chist = np.empty((t_res, ndim, ndim))
	CLV = np.empty((t_N-t_0+1, ndim, ndim))
	Mhist = np.empty((t_N-t_0+1, ndim, ndim))

#--------------------------------------------------------------------------------
#							 Forward Integration
#						Rhist[i,:,:] is growth from i to i+1
#						Qhist[i,:,:] is at i+1
#--------------------------------------------------------------------------------
	for i in range(tmax):
		M = np.dot(tlm[i,:,:],M)
		if (np.mod(i+1,rescale)==0):
			Q,R = LA.qr(M)
			Qhist[cnt,:,:] = Q
			Rhist[cnt,:,:] = R
			det=LA.det(R)
#			print("Determinant is",i, det)
			R_inv = LA.inv(R)
			M = Q
			rdiag = abs(np.diag(R))
			lyap_sum =lyap_sum+np.log(rdiag)/(rescale*dt/oneday)
			cnt = cnt+1

#--------------------------------------------------------------------------------
#							Computing Lyapunov Exponents	 
#--------------------------------------------------------------------------------
	LE = lyap_sum/cnt
#	print (LE)
#	np.savetxt('LE_takuma_ic.txt',LE)
	plt.plot(LE)
	plt.savefig("LE_takuma_ic.png",dpi = 300)
	print("Lyapunov Exponents Calculated ...")
# Initializing C matrix
	C = np.random.normal(0.0,1.0,(ndim,ndim))
	C = normalize_column(np.triu(C))

#--------------------------------------------------------------------------------
#								Backward Integration
#--------------------------------------------------------------------------------
#	print(np.any(np.isnan(C)))
	isprint =0
	for i in reversed(np.arange(t_0,t_res)):
#		Chist[i,:,:]=normalize(C,norm='l2',axis=0)
		Chist[i,:,:]=normalize_column(C)
		if (isprint==0):
			print(i)
			isprint =1
		R=Rhist[i,:,:]
#		print("Determinants is",LA.det(R))
		R_inv = LA.inv(R)
		C = np.dot(R_inv,C)
	print("Backward Integration done")	

#--------------------------------------------------------------------------------
#							Calculating Lyapunov Vectors
#--------------------------------------------------------------------------------
	for idx,i in enumerate(np.arange(t_0,t_N)):
		Q = Qhist[i-1,:,:]
		C = Chist[i,:,:]
		CLV[idx,:,:] = np.dot(Q,C)
		Mhist[idx,:,:] = tlm[i,:,:]
	print("Calculating Lyapunov Vector done")	


#--------------------------------------------------------------------------------
#								Loading/Saving data
#--------------------------------------------------------------------------------
#	np.save("CLV_takuma_ic.npy",CLV)
#	np.save("Mhist_takuma_ic.npy",Mhist)
#	CLV=np.load("CLV.npy")
#	Mhist=np.load("Mhist.npy")
#	print("Loaded")

#--------------------------------------------------------------------------------
#			 Calculating growth rates of the error along different CLVs
#--------------------------------------------------------------------------------
	CLV_list=[0,1,4,9,14,19,24,29] #CLV number to compute growth rate for
	del_t = 1000 # time to compute growth rates for
	norm_vi = np.empty((t_N-t_0-del_t+1,del_t+1,len(CLV_list)))
	
	for idx,CLV_num in enumerate(CLV_list):
		for t_index in range(t_N-t_0+1-del_t):
			vi = CLV[t_index,:,CLV_num]
			assert LA.norm(vi)>0.0
			norm_vi[t_index,0,idx] = LA.norm(vi)
			for dtx in range(del_t):
				vi = np.dot(Mhist[t_index+dtx,:,:],vi[:])
				assert LA.norm(vi)>0.0
				norm_vi[t_index,dtx+1,idx] = LA.norm(vi)
	
	norm_log_mean = np.mean(np.log(norm_vi),axis = 0)
	norm_mean_log = np.log(np.mean(norm_vi,axis=0))
	growth_log_mean = np.copy(norm_log_mean)
	growth_mean_log = np.copy(norm_mean_log)
	
	for dtx in reversed(range(del_t)):
		growth_log_mean[dtx+1,:] -= growth_log_mean[dtx,:]
		growth_mean_log[dtx+1,:] -= growth_mean_log[dtx,:]
	print("Calculating error Growth... done")	
	

#--------------------------------------------------------------------------------
#			Plotting growth rates of the error along different CLVs
#--------------------------------------------------------------------------------
	x = (np.arange(float(del_t)+0.5))*(dt/oneday)
	print(np.shape(x))
	print(np.shape(growth_log_mean))
	for idx,CLV_num in enumerate(CLV_list):
		plt.plot(x,growth_log_mean[:,idx], label="CLV num%d" %(CLV_num+1))
	plt.xlabel("days")
	plt.legend()
	plt.savefig("fig8_growth_log_mean.png",dpi=300)
	plt.close()


	for idx,CLV_num in enumerate(CLV_list):
		plt.plot(x,growth_mean_log[:,idx], label="CLV num%d" %(CLV_num+1))
	plt.xlabel("days")
	plt.legend()
	plt.savefig("fig8_growth_mean_log.png",dpi=300)
	plt.close()


#--------------------------------------------------------------------------------
#							Defining other functions	
#--------------------------------------------------------------------------------
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
