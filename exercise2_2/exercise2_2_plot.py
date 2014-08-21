from numpy import *
import matplotlib.pyplot as plt

i,diff = loadtxt('difference.dat', unpack=True)
plt.figure()
plt.loglog(10**i,diff,'ro', 10**i,diff,'-b')
plt.xlabel('i')
plt.ylabel('difference')
plt.show()