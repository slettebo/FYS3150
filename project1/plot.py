from numpy import *
import matplotlib.pyplot as plt

v = loadtxt('/home/mammothopus/Desktop/FYS3150 - Computational Physics/GIT/FYS3150/build-project1-Desktop-Debug/u_numerical_solution_10.dat', unpack=True)
n = len(v)
h = 1./(n+1)
x = linspace(0,1,n)
u = 1 - (1-exp(-10))*x - exp(-10*x)
f = 100*exp(-10*x)*h*h;

print n
plt.figure()
plt.plot(x,u,'-g',x,v,'-r')
plt.legend(['Analytical', 'Numerical'])
plt.xlabel('x')
plt.ylabel('y')
plt.show()
