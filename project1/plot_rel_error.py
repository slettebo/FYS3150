from numpy import *
import matplotlib.pyplot as plt


n = array([10,100,1000,10000,100000,1000000])
rel_error = zeros(len(n))

for i in range(0,len(n)):
	rel_max = (loadtxt('/home/mammothopus/Desktop/FYS3150 - Computational Physics/GIT/FYS3150/build-project1-Desktop-Debug/relative_error_%d.dat' %n[i], unpack=True)).max()
	rel_error[i] = rel_max

h = 1./(n+1)

print h
print rel_error

plt.figure()
plt.grid('on')
plt.title('Plot of $log_{10}(| \\frac{v_i - u_i}{u_i} |)$ vs. $log_{10}(h)$: rel.error vs. step length')
plt.plot(log10(h),rel_error,'--k',log10(h),rel_error,'or')
#plt.legend(['Analytical', 'Numerical'])
plt.xlabel('$log_{10}(h)$')
plt.ylabel('$log10(\\epsilon)$')
plt.show()
