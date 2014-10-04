
# Exercise 2

import stagger
import matplotlib.pylab as plt
from numpy import *

N = 1024
T = 10
dt = 0.00001

z   = linspace(0,1,N)
u   = zeros((N,T))
rho = zeros((N,T))

rho[where((z < 0.6) & (z > 0.3)), 0] = 1.0

plt.plot(z,rho[:,0])
plt.show()

dz  = z[1]-z[0]

for i in xrange(T-1):
    rho[:,i+1] = rho[:,i] \
            - dt * (u[:,i] * stagger.deriv( stagger.interp(rho[:,i],direction=-1), dz) \
            + rho[:,i] * stagger.deriv( stagger.interp(u[:,i],direction=-1), dz) )

    u[:, i+1] = u[:,i] \
            - stagger.interp(u[:,i]) * stagger.deriv(u[:,i], dz) \
            - 5./3 * rho[:,i] ** 2./3 \
            * stagger.deriv( stagger.interp( rho[:,i], direction=-1), dz)


plt.plot(z,rho[:,2])
plt.show()


plt.plot(z,rho[:,5])
plt.show()
