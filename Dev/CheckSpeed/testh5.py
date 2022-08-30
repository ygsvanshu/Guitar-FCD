import numpy
import h5py

f = h5py.File('test.h5','r')
print(f['M'][:])