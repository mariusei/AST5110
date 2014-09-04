# AST5110
# Thu 04/09/2014
# Marius Berge Eide

import matplotlib.pylab as plt
import matplotlib.animation as animation
from numpy import *

# Solve advection eqs

# Initialise

class advection:
    def __init__(self, N = 100, vz  = 1.0, mname = 'FTCS'):
         
        self.N  = N
        self.vz = vz 
        self.mname = mname 

        z = linspace(0,1,N, endpoint=False)   # not inclusive endpoint

        rho     = zeros(N)
        rho_next= zeros(N)

        condition = where((z > 0.2) & (z < 0.3))

        rho[:]  = 0.1
        rho[condition] = 1.0

        # FTCS to solve advection eq
        # Forward in time
        # Centered in space

        t = 0
        dt = 0.001        # random
        dz = z[1]-z[0]  # assuming constant z step length
        delta = vz * dt / (2.0 * dz)

        # boundary conditions
        # Assume locked boundaries
        rho_next[0]     = rho[0]
        rho_next[N-1]   = rho[N-1]

        # make these available to all class
        self.z          = z
        self.rho        = rho
        self.rho_next   = rho_next
        self.delta      = delta

        if mname == 'staggered_leapfrog'

    def step(self):
        rho_next, rho = self.rho_next, self.rho
        
        if self.mname   == 'FTCS':
            for iz in xrange(1,self.N-1):
                self.rho_next[iz] = rho[iz] \
                                    - self.delta * ( rho[iz+1] - rho[iz-1] )
        elif self.mname == 'lax':
            for iz in xrange(1,self.N-1):
                self.rho_next[iz] = 0.5 * (rho[iz+1] + rho[iz-1]) \
                        - self.delta * ( rho[iz+1] - rho[iz-1] )

        elif self.mname == 'staggered_leapfrog':
            

        else:
            raise NameError('Wrong method specified: '+str(self.mname))

        self.rho = rho_next

    def get_state(self):
        return (self.z, self.rho)


ftcs = advection(mname='FTCS')
lax  = advection(mname='lax')  # independent, such as not to mix rho-arrays

# Animation

for method in [ftcs, lax]:
    fig = plt.figure()
    minval = min(method.get_state()[0])
    maxval = max(method.get_state()[1])

    ax  = plt.axes(xlim=(minval, maxval), ylim=(-0.5,2))
    line,   = ax.plot([], [], lw=2)


    def animate(it):
        method.step()
        line.set_data(method.get_state())
        return line,

    def init():
        line.set_data([], [])
        return line,

    anim = animation.FuncAnimation(fig, animate, init_func=init, \
            frames=300, interval=20, blit=True)

    anim.save(method.mname+'.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    plt.show()
