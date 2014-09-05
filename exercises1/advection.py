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
        
        # Initialising, using self.variable where self makes the variable
        # available to the entire instance of the class
        # e.g. you initialise an instance: 
        #   mynameofinstance = ADVECTION(mname='staggered_leapfrog')
        # then you can retrieve N as
        #   mynameofinstance.N
        # or call functions
        #   mynameofinstance.step()
        
        self.N  = N
        self.vz = vz 
        self.mname = mname 

        z = linspace(0,1,N, endpoint=False)   # not inclusive endpoint

        rho     = zeros(N)
        rho_next= zeros(N)

        condition = where((z > 0.2) & (z < 0.3))

        rho[:]  = 0.1
        rho[condition] = 1.0

        t = 0
        dt = 0.001      # random
        dz = z[1]-z[0]  # assuming constant z step length
        delta = vz * dt / (2.0 * dz)

        # boundary conditions
        # Assume locked boundaries
        rho_next[0]     = rho[0]
        rho_next[N-1]   = rho[N-1]

        # make these available to entire class
        self.z          = z
        self.rho        = rho
        self.rho_next   = rho_next
        self.delta      = delta

        if mname == 'staggered_leapfrog':
            # Construct N * 2 dim array to hold time step
            # n-2 and n-1 needed to find time step n
            self.rho_tmp = zeros((N,2), float)
            self.rho_tmp[:,0] = rho
            self.rho_tmp[:,1] = rho
            self.delta = vz * dt / dz
        if mname == 'upwind':
            self.delta = vz * dt / dz

    def step(self):
        rho_next, rho = self.rho_next, self.rho
        if self.mname == 'staggered_leapfrog':
            rho_tmp = self.rho_tmp
        
        if self.mname   == 'FTCS':
            for iz in xrange(1,self.N-1):
                self.rho_next[iz] = rho[iz] \
                                    - self.delta * ( rho[iz+1] - rho[iz-1] )
        elif self.mname == 'lax':
            for iz in xrange(1,self.N-1):
                self.rho_next[iz] = 0.5 * (rho[iz+1] + rho[iz-1]) \
                        - self.delta * ( rho[iz+1] - rho[iz-1] )

        elif self.mname == 'staggered_leapfrog':
            for iz in xrange(1,self.N-1):
                self.rho_next[iz] = rho_tmp[iz,0] \
                        - self.delta * ( rho_tmp[iz+1,1] - rho_tmp[iz-1,1])
            self.rho_tmp[:,0] = rho_tmp[:,1]
            self.rho_tmp[:,1] = rho_next

        elif self.mname == 'upwind':
            for iz in xrange(1,self.N-1):
                if self.vz > 0:
                    self.rho_next[iz] = rho[iz] \
                            - self.delta * ( rho[iz] - rho[iz-1] )
                else:
                    self.rho_next[iz] = rho[iz] \
                            - self.delta * ( rho[iz+1] - rho[iz] )

        else:
            raise NameError('Wrong method specified: '+str(self.mname))

        self.rho = rho_next

    def get_state(self):
        return (self.z, self.rho)

# initialise instances of class where mname gives the choice of model
ftcs = advection(mname='FTCS')
lax  = advection(mname='lax')  # independent, such as not to mix rho-arrays
s_leap=advection(mname='staggered_leapfrog')
upwind=advection(mname='upwind')

# Animation
# iterate through methods and produce animations

for method in [ftcs, lax, s_leap, upwind]:
#for method in [upwind]:
#for method in [s_leap]:
    fig = plt.figure()
    minval = min(method.get_state()[0])
    maxval = max(method.get_state()[1])

    # set up axes and plot empty data set
    ax  = plt.axes(xlim=(minval, maxval), ylim=(-0.5,2))
    line,   = ax.plot([], [], lw=2)


    def animate(it):
        method.step()
        line.set_data(method.get_state())
        return line,

    def init():
        line.set_data([], [])
        return line,

    # Do the real animation part
    # where animation comes from matplotlib.animation
    # and FuncAnmation requires figure fig, functions animate and init
    # and blit means only redraw parts that are changing
    anim = animation.FuncAnimation(fig, animate, init_func=init, \
            frames=300, interval=20, blit=True)

    # To save, uncomment following line
    # Requires packages ffmpeg, mencoder and libx264 with dependencies
    #anim.save(method.mname+'.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    # the show goes on forever:
    plt.show()
