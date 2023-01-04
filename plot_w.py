import numpy as np
import sys
import matplotlib.pyplot as plt
import h5py

N=sys.argv[1]
sys_size=int(sys.argv[2])

f=h5py.File("field.h5",'r')

w2=f.get(N).value
w=np.reshape(w2,(sys_size,sys_size+2))
w=np.delete(w,[sys_size,sys_size+1],1)

plt.imshow(w)
plt.show()
