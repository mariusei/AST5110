
from numpy import *
import matplotlib.pylab as plt

def forwardsubst(matrix, b):
    
    M       = max(shape(matrix))
    xout    = zeros(M)

    #for i in xrange(0,M):
    #    xout[i] = ( b[i] - dot(matrix[i, 0:i], xinit[0:i]) ) / matrix[i, i]

    xout[0] = 1./matrix[0,1] * b[0]
    #for i in xrange(1,M):
    #    xout[i] = ( b[i] - dot(matrix[i], b[i-1:i+1]) ) / matrix[i,1]

    for i in xrange(1,M):
        xout[i] = ##########


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

N   = 512
x   = linspace(0,1,N)
rhs = zeros(N)
phi0= zeros(N)
phi1= zeros(N)

# RHS
h       = x[1]-x[0]
b       = h**2/2. * sourceterm(x)
uvec    = ones(N-1)     # upper diagonal of 1's
phi0[:] = 1.            # initial guess of solution

rhs = b 
rhs[0:N-1] += uvec * phi0[0:N-1]


# LHS
lvec    = zeros((N,2))  # lower diagonal
lvec[1:,0]=  1
lvec[:,1] = -2


# Iterate
for j in xrange(2):
    phi1[:] = forwardsubst(lvec, rhs)
    rhs[:]  = b
    rhs[0:N-1]+= uvec * phi1[0:N-1]

plt.plot(x, phi1)
plt.show()

