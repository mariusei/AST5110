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


def eos(rho, gamma):
    return power(rho, gamma)

def ecsolver(f, dx, dt, xi, yi):
    """ Advances one step in time """
    y1 = yi + dt * f(xi, yi, dx, dt)
    
    rho  = yi[0]
    p    = yi[1]
    p_g  = yi[2]
    Q    = yi[3]
    #Q    = deriv( deriv( deriv( deriv(p, dx, direction=-1), dx), dx, direction=-1), dx)

    #y1[4] =     -       p / interp(p, direction=1) \
    #                   *deriv( interp(p, direction=-1), dx)  \
    y1[4] =       -       deriv((p_g + Q), dx) \
    #            +       interp(rho, direction=1) * g_z

    return y1

def rk4solver(f, dx, dt, xi, yi):
    """ Advances one step in time using a predictor-corrector method """

    k1 = f(xi, yi, dx, dt)
    k2 = f(xi, yi + dt/2 * k1, dx, dt/2)
    k3 = f(xi, yi + dt/2 * k2, dx, dt/2)
    k4 = f(xi, yi + dt   * k3, dx, dt)

    y1 = yi + dt/6 * (k1 + 2.*k2 + 2.*k3 + k4)
    
    rho  = yi[0]
    p    = yi[1]
    p_g  = yi[2]
    Q    = yi[3]
    #Q    = deriv( deriv( deriv( deriv(p, dx, direction=-1), dx), dx, direction=-1), dx)

    y1[4] =     -       p / interp(p, direction=1) \
                       *deriv( interp(p, direction=-1), dx)  \
           -       deriv((p_g + Q), dx) \
                +       interp(rho, direction=1) * g_z

    return y1

def cont_eqs(xi, yi, dx, dt):

    rho  = yi[0]
    p    = yi[1]
    p_g  = yi[2]
    Q    = yi[3]

    gamma = 5./3
    
    p_g  = 0 #eos(rho, gamma)

    dpdt =      -       p / interp(p, direction=1) \
                    *deriv( interp(p, direction=-1), dx) \
        -       deriv((p_g + Q), dx) \
                +       interp(rho, direction=1) * g_z

    drhodt =    -       deriv( p, dx)

    dp_gdt =    zeros(N) #          gamma * power(rho, gamma-1) * drhodt

    return array([drhodt, dpdt, dp_gdt, zeros(N),    zeros(N)])

### ANIMATION ###


def get_state(yinout):

    yinout[:] = rk4solver(f=cont_eqs, dx=dx, dt=dt, xi=0, \
            yi=yinout)

    return (x, yinout)

# y is THE array, holding values that will be integrated in time

N = 1024
dt= 0.001

x = linspace(0,1,N)
dx= x[1] - x[0]

rho  = ones(N)
#rhov[where((x < 0.7) & (x > 0.4))] = 2.0
rhov = ones(N) + peak(N/2, N/6, N, 3)
p_g  = zeros(N) #eos(rho, 5./3)
Q    = zeros(N)
g_z  = 0.0

y    = zeros((5,N)) ####### 
y[0] = rho
y[1] = rhov
y[2] = p_g
y[3] = Q


fig = plt.figure()
minval = min(x)
maxval = max(x)

# set up axes and plot empty data set
ax  = plt.axes(xlim=(minval, maxval), ylim=(-2.5,13))
line_rho,   = ax.plot([], [], lw=2)
line_rv,    = ax.plot([], [], lw=2)
line_p_g,   = ax.plot([], [], lw=2)
line_diffop,   = ax.plot([], [], lw=2)
lines = [line_rho, line_rv, line_p_g, line_diffop]

def init():
    for line in lines:
        line.set_data([], [])
    ax.legend(lines, \
            [r'$\rho$', r'$\rho v$', r'$p_{\rm gas}$', 'diffoperator'])
    return lines

def animate(it, y, lines):
    x, y = get_state(y)
    rho, rhov, p_g, Q   = y[0], y[1], y[2], y[3]
    diffterm = y[4] ######
    lines[0].set_data(x, rho)
    lines[1].set_data(x, rhov)
    lines[2].set_data(x, p_g)
    lines[3].set_data(x, diffterm)
    return lines


anim = animation.FuncAnimation(fig, animate, init_func=init, \
        fargs=(y, lines), frames=1024, interval=05, blit=True)

# To save, uncomment following line
# Requires packages ffmpeg, mencoder and libx264 with dependencies
print 'Writing animation to file...'
#anim.save('no_press_rhov.mp4', fps=125, extra_args=['-vcodec', 'libx264'])
print 'Done'

# the show goes on forever:
plt.show()







## TEST SUITE:
#
## Test time integrator

#def linear(xi, yi, dx, dt):
#    return yi
#
#N = 1024
#tmparr = zeros(N)
#ytmp =  zeros((N,2))
#ytmp[0,:] = 1.0
#for i in xrange(N-1):
#    ytmp[i+1,0] = ecsolver(linear, dx,  10*dt, 0.0, ytmp[i,0])
#    ytmp[i+1,1] = rk4solver(linear, dx, 10*dt, 0.0, ytmp[i,1])
#    tmparr[i+1] = tmparr[i] + 10*dt
#
#plt.plot(tmparr, ytmp[:,0])
#plt.plot(tmparr, ytmp[:,1])
#plt.plot(tmparr, exp(tmparr))
#plt.show()
#
#
