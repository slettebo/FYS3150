from numpy import *
import matplotlib.pyplot as plt

for i in [10,100,1000]:
    v = loadtxt('/home/mammothopus/Desktop/FYS3150 - Computational Physics/GIT/FYS3150/build-project1-Desktop-Debug/u_numerical_solution_%d.dat' %i, unpack=True)
    n = len(v)
    h = 1./(n+1)
    x = linspace(0,1,n)
    u = 1 - (1-exp(-10))*x - exp(-10*x)
    f = 100*exp(-10*x)*h*h;

    print n
    plt.figure()
    plt.title('n=%d' %i)
    plt.plot(x,u,'--b',x,v,'-r')
    plt.legend(['Analytical', 'Numerical'])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
