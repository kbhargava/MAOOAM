#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import linalg as LA
ndim =36

oneday =8.64

def main():
    LE=read_file()
    plt.plot(LE)
    plt.savefig("trial.png",dpi=300)
def read_file():
    LE=np.arange(36)
    return LE

if __name__ == "__main__":
    main()
