#!/usr/bin/env python

# Exercise 2

import stagger
import matplotlib.pylab as plt
import matplotlib.animation as animation
from numpy import *


def eos(rho, gamma):
    return power(rho, gamma)

class MomEq:

    def __init__(self):
        self.mname = 'Momentum Eq. Stagger solved'

        N = 1024
        dt = 0.0001

        z   = linspace(0,1,N)
        rho = ones((N,2))
        rhov= zeros((N,2))
        p_g = zeros(N)
        Q   = zeros(N)

        rho[:,0] = 1.0 #1.0 + self.peak(z, N)
        rhov[:,0] = 0  + self.peak(z,N) # * stagger.interp(rho[:,0], direction=-1)
 #       rhov[where((z>0.3) & (z< 0.5)),0] = 1.0

        dz  = z[0] - z[1]
        g_z = 0.0

        self.z, self.rho, self.rhov = z, rho, rhov
        self.Q, self.p_g    = Q, p_g
        self.dt, self.dz    = dt, dz
        self.g_z            = g_z

        self.a  = zeros(N)
        self.b  = zeros(N)

    def peak(self, x, N):

        y = zeros(N)
        y[N/2 - 2*N/16 : N/2 + 2*N/16] = 1.0
        yy = 0.002*sin(x[0:N/4] * 4 * pi)

        return convolve(y, yy, 'same')

    def eos(self, gamma):
        self.p_g = power(self.rho[:,0], gamma)

    def step(self, it):
        rhov, rho, Q    = self.rhov, self.rho, self.Q
        dt, dz, g_z     = self.dt, self.dz, self.g_z

        self.eos(5./3)
        p_g = self.p_g

        a   =   - stagger.deriv(stagger.interp(rhov[:,0]**2, direction=-1)\
                /rho[:,0], dz) + roll(rhov[:,0],+300)
        b   =   0 #- stagger.deriv(( \
                #        p_g + Q \
                #        ), dz)

        c   =   stagger.interp(rho[:,0], direction=1) * g_z 
        
        rhov[:,1] = rhov[:,0] + dt \
                 * (a + b + c)
#                * (- 1./(stagger.interp(rho[:,0], direction=1)) \
#                    * rhov[:,0] \
#                    * stagger.deriv(stagger.interp(rhov[:,0], direction=-1), dz) \
#                 - stagger.deriv(( \
#                        p_g + Q \
#                        ), dz) \
#                 + stagger.interp(rho[:,0], direction=1) * g_z )

        self.a, self.b, self.c = a,b,c

        rho[:,1] = rho[:,0] -  dt * stagger.deriv( rhov[:,1], dz, direction=-1)

        self.rhov[:,0] = rhov[:,1]
        self.rho[:,0] = rho[:,1]
        

    def get_state(self):
        return (self.z, self.a) # self.rhov[:,0])
        
        


method = MomEq()

method


fig = plt.figure()
minval = min(method.get_state()[0])
maxval = max(method.get_state()[0])

# set up axes and plot empty data set
ax  = plt.axes(xlim=(minval, maxval), ylim=(-0.5,2))
ax.set_title('Method: ' + method.mname)
line,   = ax.plot([], [], lw=2)


def animate(it):
    method.step(it)
    line.set_data(method.get_state())
    ax.set_title('it=%g'%it)
    return line,

def init():
    line.set_data([], [])
    return line,

# Do the real animation part
# where animation comes from matplotlib.animation
# and FuncAnmation requires figure fig, functions animate and init
# and blit means only redraw parts that are changing
anim = animation.FuncAnimation(fig, animate, init_func=init, \
        frames=60, interval=10, blit=True)

# To save, uncomment following line
# Requires packages ffmpeg, mencoder and libx264 with dependencies
print 'Writing animation to file...'
#anim.save(method.mname+'.mp4', fps=15, extra_args=['-vcodec', 'libx264'])
print 'Done'

# the show goes on forever:
plt.show()
