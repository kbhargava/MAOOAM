#!/usr/bin/python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import linalg as LA
ndim =36

oneday =8.64

def main():
    traj, tlm, tim = read_file("evol_field_tlm_takuma_ic_2.dat")
    print("Loaded")
    dt = tim[1]-tim[0]
    print(dt)
    tmax = len(tim)
    rescale = 10
    
    lyap_sum =np.zeros(ndim)
    cnt = 0
    M = np.eye(ndim)
    S = np.eye(ndim)
    LE = np.empty(ndim)
    
    for i in range(0,tmax):
        M = np.dot((tlm[i,:,:]),M)
        if (np.mod(i+1,rescale)==0):
            Q,R = LA.qr(M)
            M = Q
            rdiag=abs(np.diag(R))
            lyap_sum =lyap_sum+np.log(rdiag)/(rescale*dt/oneday)
            cnt = cnt+1
    LE = lyap_sum/cnt
#    S = np.dot(M.T,M)
#    w,v = LA.eig(S)
#    LE=np.log(abs(w))/(2*tmax/oneday)
#    LE=np.sort(LE)
#    print (LE)
#    LE=np.flip(LE,axis=0)
    print (LE)
#    np.savetxt('LE_1200.txt',LE)
    plt.plot(LE)
    plt.savefig("LE.png",dpi = 300)



def read_file(ifile):
    #return np.ndarray[time, ndim]
    with open(ifile,"r") as f:
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
