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
dt= 0.0001

x = linspace(0,1,N)
dx= x[1] - x[0]

rho  = ones(N)
rhov = ones(N) + peak(N/2, N/6, N, 3)
p_g  = zeros(N)
Q    = zeros(N)
g_z  = 0.0


rr = zeros(N)
rv = zeros(N)

def step_rhov_rho(dx, y):

    rhov, rho = y[0], y[1]

    p_g[:] =  0 #0.01* power(rho, 5./3)

    return array([-    rhov[:] / interp(rho[:], direction=1) \
                    *deriv( interp(rhov[:]), dx) \
                -    deriv( p_g + Q, dx) \
                +    interp(rho[:], direction=1) * g_z \
            ,    - deriv(rhov[:], dx) ])

def step(t):

#    p_g[:] = zeros(N) #power(rho, 5./3)
#
#    rv[:] = rhov[:] + dt * ( \
#                -    rhov[:] / interp(rho[:], direction=1) \
#                    *deriv( interp(rhov[:]), dx) \
#                -    deriv( p_g + Q, dx) \
#                +    interp(rho[:], direction=1) * g_z )
#    
#    rr[:] = rho[:] + dt * ( - deriv(rhov[:], dx) )

#    rv1, rr1   = step_rhov_rho(dx,     [rhov, rho])
#    rv2, rr2   = step_rhov_rho(0.5*dx, [rhov +0.5*rv1, rho+0.5*rr1])
#    rv3, rr3   = step_rhov_rho(0.5*dx, [rhov +0.5*rv2, rho+0.5*rr2])
#    rv4, rr4   = step_rhov_rho(dx,     [rhov +    rv3, rho+    rr3])

    rv, rr = step_rhov_rho(dx, [rhov, rho])

    rhov[:] = rhov[:] + dt * rv # + dt/6 * (rv1 + 2.*rv2 + 2.*rv3 + rv4)
    rho[:]  = rho[:]  + dt * rr # + dt/6 * (rr1 + 2.*rr2 + 2.*rr3 + rr4) 

#    rhov[:] = rv[:]
#    rho[:]  = rr[:]



### ANIMATION ###

def get_state():
    return (x, rhov, x, rho, p_g)

def animate(it):
    step(it)
    line_rv.set_data(get_state()[0], get_state()[1])
    line_rho.set_data(get_state()[2], get_state()[3])
    line_p_g.set_data(get_state()[0], get_state()[4])
    ax.set_title('it=%g'%it)
    return line_rv, line_rho, line_p_g

def init():
    line_rho.set_data([], [])
    line_rv.set_data([], [])
    line_p_g.set_data([], [])
    ax.legend([line_rho, line_rv, line_p_g], [r'$\rho$', r'$\rho v$', r'$p_{\rm gas}$'])
    return line_rho, line_rv, line_p_g

fig = plt.figure()
minval = min(get_state()[0])
maxval = max(get_state()[0])

# set up axes and plot empty data set
ax  = plt.axes(xlim=(minval, maxval), ylim=(-0.5,5))
line_rho,   = ax.plot([], [], lw=2)
line_rv,    = ax.plot([], [], lw=2)
line_p_g,   = ax.plot([], [], lw=2)


anim = animation.FuncAnimation(fig, animate, init_func=init, \
        frames=10240, interval=15, blit=True)

# To save, uncomment following line
# Requires packages ffmpeg, mencoder and libx264 with dependencies
print 'Writing animation to file...'
#anim.save('alive.mp4', fps=480, extra_args=['-vcodec', 'libx264'])
print 'Done'

# the show goes on forever:
plt.show()
