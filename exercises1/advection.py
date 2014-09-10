# AST5110
# Thu 04/09/2014
# Marius Berge Eide

import matplotlib.pylab as plt
import matplotlib.animation as animation
from numpy import *

# Solve advection eqs

# Initialise

class ADVECTION: 
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

        scale   = 1.0  # Scaling factor
        dz = z[1]-z[0]  # assuming constant z step length

        # boundary conditions
        # Assume locked boundaries
        rho_next[0]     = rho[0]
        rho_next[N-1]   = rho[N-1]

        # make these available to entire class
        self.z          = z
        self.rho        = rho
        self.rho_next   = rho_next

        if mname == 'FTCS':
            self.delta = scale * dz

        elif mname == 'lax':
            # Courant criterion
            dt =  dz / abs(vz)
            self.delta = vz * dt / (2.0 * dz)
            print vz * dt / dz 

        elif mname == 'staggered_leapfrog':
            # Construct N * 2 dim array to hold time step
            # n-2 and n-1 needed to find time step n
            self.rho_tmp = zeros((N,2), float)
            self.rho_tmp[:,0] = rho
            self.rho_tmp[:,1] = rho
            dt = scale * dz / abs(vz)
            self.delta = vz * dt / dz
        elif mname == 'upwind':
            dt         = scale * dz / abs(vz)
            self.delta = vz * dt / dz

        elif mname == 'van_leer':
            dt          = 0.1 * dz
            self.delta  = vz * dt / dz
            # Construct array with twice as many points
            self.z   = linspace(0,1,2*N, endpoint=False)
            self.rho = zeros((2*N), float)
            self.rho[:] = 0.1
            self.rho[where((self.z > 0.2) & (self.z < 0.3))] = 1.0
            self.rho_next = self.rho

        else:
            raise ValueError('Unknown method: %s.'%mname)



    def step(self, it):
        rho_next, rho = self.rho_next, self.rho
        delta         = self.delta
        if self.mname == 'staggered_leapfrog':
            rho_tmp = self.rho_tmp
        
        if self.mname   == 'FTCS':
            for iz in xrange(1,self.N-1):
                rho_next[iz] = rho[iz] \
                                    - delta * ( rho[iz+1] - rho[iz-1] )
        elif self.mname == 'lax':
            for iz in xrange(1,self.N-1):
                rho_next[iz] = 0.5 * (rho[iz+1] + rho[iz-1]) \
                        - delta * ( rho[iz+1] - rho[iz-1] )
            if it == 2:
                print self.rho
                print rho_next
                print rho, delta,self.N

        elif self.mname == 'staggered_leapfrog':
            for iz in xrange(1,self.N-1):
                rho_next[iz] = rho_tmp[iz,0] \
                        - delta * ( rho_tmp[iz+1,1] - rho_tmp[iz-1,1])
            self.rho_tmp[:,0] = rho_tmp[:,1]
            self.rho_tmp[:,1] = rho_next

        elif self.mname == 'upwind':
            for iz in xrange(1,self.N-1):
                if self.vz > 0:
                    rho_next[iz] = rho[iz] \
                            - delta * ( rho[iz] - rho[iz-1] )
                else:
                    rho_next[iz] = rho[iz] \
                            - delta * ( rho[iz+1] - rho[iz] )

        elif self.mname == 'van_leer':
            for iz in xrange(2,self.N*2-3):
                rho_next[iz] = rho[iz] \
                        - delta * (rho[iz+1] - rho[iz-1] \
                        + 0.5*(1-delta) * (self.slope(rho,iz+1)*rho[iz] \
                                          -self.slope(rho,iz-1)*rho[iz]))
                
            

        else:
            raise NameError('Wrong method specified: '+str(self.mname))

        self.rho = rho_next

    def get_state(self):
        return (self.z, self.rho)

    def slope(self,rho,i):
        if ( (rho[i] - rho[i-2]) * (rho[i+2] - rho[i]) > 0 ):
            return 2 * (rho[i]   - rho[i-2]) * (rho[i+2] - rho[i]) \
                     / (rho[i+2] - rho[i-2])
        else:
            return 0

# initialise instances of class where mname gives the choice of model
ftcs        = ADVECTION(mname='FTCS')
lax         = ADVECTION(vz=1.05, mname='lax')  
s_leap      = ADVECTION(mname='staggered_leapfrog')
upwind      = ADVECTION(mname='upwind')
van_leer    = ADVECTION(mname='van_leer')

# Animation
# iterate through methods and produce animations

#for method in [ftcs, lax, s_leap, upwind]:
#for method in [upwind]:
for method in [lax]:
#for method in [s_leap]:
#for method in [van_leer]:
    fig = plt.figure()
    minval = min(method.get_state()[0])
    maxval = max(method.get_state()[1])

    # set up axes and plot empty data set
    ax  = plt.axes(xlim=(minval, maxval), ylim=(-0.5,2))
    line,   = ax.plot([], [], lw=2)


    def animate(it):
        method.step(it)
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

    print method.delta
