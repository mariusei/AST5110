
from numpy import *
import matplotlib.pylab as plt
import stagger


# Exercise set 3:
# Why Multigrid Is Needed
# From Claudio
# AST5110 H2014
# Marius Berge Eide
#
# Solve:
# d2 phi / dx2 = rho(x)
# 
# 1) Using Gauss-Seidel, lexiographic ordering
# 2) Find residuals (d2 phi / dx2 - rho(x) = 0)
# 3) Find error (solve eq analytically with source function
#   rho(x) = ... in sourceterm
#   with solution in solutionterm
# 4) Understand that the choice of grid size determines what eigenvectors
#    will dominate the solver
#    where smaller grid sizes tend to favour smaller oscillations
#    and lager grid sizes tend to favour larger structures
#

def forwardsubst(upper, lower, diagonal, bvec, xin):
    
    M       = max(shape(diagonal))
    xout    = zeros(M)

    # PERIODIC BOUNDARY CONDITIONS: ELEMENTS 0M = 1 and M0 = 1
    xout[0] = 1./diagonal[0] * ( bvec[0] - upper[0] * xin[1] \
                                +1.0 * xin[M-1] )
    for i in xrange(1,M-1):
        xout[i] = 1./diagonal[i] * ( bvec[i]   - lower[i-1] * xout[i-1] \
                                            - upper[i] * xin[i+1] )
    xout[M-1] = 1./diagonal[M-1] * ( bvec[M-1] - lower[M-2] * xout[M-1]    \
                                  +1.0 * xout[0] )


    return xout

def sourceterm(x):

    a1 = 500
    a2 = 50
    k1 = 20
    k2 = 1

    return a1 * sin(2*pi* k1* x)  +  a2 * sin(2*pi* k2* x) 

def solutionterm(x):

    a1 = 500
    a2 = 50
    k1 = 20
    k2 = 1

    return  -  a1 / (2*pi*k1)**2 * sin(2*pi* k1* x) \
            -  a2 / (2*pi*k2)**2 * sin(2*pi* k2* x) 

# Number of nodes
N   = [16, 128, 512]


fig1, axes = plt.subplots(len(N), 1)
fig2, ax2   = plt.subplots(1, 1)
fig3, ax3   = plt.subplots(1, 1)

for k in xrange(len(N)):
    n = N[k]

    x   = linspace(0,1,n)
    phi0= zeros(n)
    phi1= zeros(n)
    solu= solutionterm(x)

    M = 1024
    resid=zeros(M)
    error=zeros(M)

    # RHS
    h       = x[1]-x[0]
    b       = h**2/2. * sourceterm(x)


    # Iterate
    for j in xrange(M):
        phi1[:] = forwardsubst(ones(n-1), ones(n-1), -2*ones(n), b, phi0[:])
        phi0[:] = phi1[:]
        error[j] = sum( abs( phi0 - solu ) ) / (n - 1)
        # Residuals: 
        # Find d2phi / dx^2 - rho(x) = 0
        # using either:
        # STAGGER 2ND DERIV:
        #resid[j] = sum( abs( stagger.deriv( stagger.deriv( \
        #                    phi1, h, direction=-1), h, direction=+1) \
        # 2ND DERIV BY HAND:
        resid[j] = sum( abs( \
                2.*( roll(phi1, -1) -2.*phi1 + roll(phi1,+1))/h**2 \
                - sourceterm(x) ) )

        # Plot four times
        if j%(M/4) == 0:
            axes[k].plot(x, phi1)
            axes[k].text(1.05, 0.5, N[k], 
                    size='large', transform=axes[k].transAxes)

    axes[k].plot(x, solu, 'k')
    ax2.plot(linspace(0,M-1,M), resid, label=('Nodes=%i' % N[k]))
    ax3.plot(linspace(0,M-1,M), error, label=('Nodes=%i' % N[k]))

ax2.set_yscale('log')
ax2.set_xscale('log')
plt.legend()
plt.show()
