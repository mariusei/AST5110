
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
        T = 10
        dt = 0.00001

        z   = linspace(0,1,N)
        rho = ones((N,2))
        rhov= zeros((N,2))
        p_g = zeros(N)
        Q   = zeros(N)

        rho[:,0] = 1.0 + self.peak(z, N)

        dz  = z[0] - z[1]
        g_z = 0.0

        self.z, self.rho, self.rhov = z, rho, rhov
        self.Q, self.p_g    = Q, p_g
        self.dt, self.dz    = dt, dz
        self.g_z            = g_z

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
        
        rhov[:,1] = rhov[:,0] + dt \
                * (- 1./rho[:,0] \
                    * stagger.interp(rhov[:,0]) \
                    * stagger.deriv(rhov[:,0], dz) \
                 - stagger.deriv(stagger.interp(( \
                        p_g + Q \
                        ), direction=-1), dz) \
                 + rho[:,0] * g_z )


        rho[:,1] = rho[:,0] -  dt * stagger.deriv( rhov[:,1], dz)
        


        self.rhov[:,0] = rho[:,1]
        self.rho[:,0] = rho[:,1]

    def get_state(self):
        return (self.z, self.rho[:,0])
        
        


method = MomEq()

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
