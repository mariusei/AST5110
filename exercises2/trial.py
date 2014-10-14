import matplotlib.pylab as plt
import matplotlib.animation as animation
from numpy import *
from stagger import deriv, interp

def peak(pos, width, N, amplitude):

    y  = zeros(N)
    yy = zeros(N)
    yy[pos-width/2:pos+width/2] = 1.0
    #peak = sin(linspace(0,pi,width))
    peak  = exp(-(linspace(0,10,width)-5)**2 / (2.*(10/8)**2))

    y = convolve(yy, peak, 'same')
    return amplitude * y / max(y)


N = 1024
dt= 0.001

x = linspace(0,1,N)
dx= x[1] - x[0]

rho  = ones(N)
rhov = ones(N) + peak(N/2, N/6, N, 3)
p_g  = zeros(N)


rr = zeros(N)
rv = zeros(N)

def step(t):

    p_g[:] = power(rho, 5./3)

    rv[:] = rhov[:] + dt * ( \
                -    rhov[:] / interp(rho[:], direction=1) \
                    *deriv( interp(rhov[:]), dz) \
                -    deriv( p_g + Q, dz) \
                +    interp(rho[:], direction=1) * g_z )
    
    rr[:] = rho[:] + dt * ( - deriv(rhov[:], dz) )

    rhov[:] = rv[:]
    rho[:]  = rr[:]

def get_state():
    return (x, rr)

def animate(it):
    step(it)
    line.set_data(get_state())
    ax.set_title('it=%g'%it)
    return line,

def init():
    line.set_data([], [])
    return line,
