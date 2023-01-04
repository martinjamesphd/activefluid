import numpy as np
import sys
import matplotlib.pyplot as plt
import h5py

N=sys.argv[1]

f=h5py.File("spectra.h5",'r')

spectrum=f.get(N).value
k=f.get("k").value

plt.loglog(k,spectrum)
plt.show()
