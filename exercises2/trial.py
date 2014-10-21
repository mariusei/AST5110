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

    return y1

def cont_eqs(xi, yi, dx, dt):

    rho  = yi[0]
    p    = yi[1]
    p_g  = yi[2]
    Q    = yi[3]

    gamma = 5./3
    
    p_g  = zeros(N) # eos(rho, gamma)

    dpdt =      -       p / interp(p, direction=1) \
                       *deriv( interp(p, direction=-1), dx) \
                -       deriv((p_g + Q), dx) \
                +       interp(rho, direction=1) * g_z

    drhodt =    -       deriv( p, dx)

    dp_gdt =            zeros(N) # gamma * power(rho, gamma-1) * drhodt

    return array([drhodt, dpdt, dp_gdt, zeros(N)])

### ANIMATION ###


def get_state(yinout):

    yinout[:] = ecsolver(f=cont_eqs, dx=dx, dt=dt, xi=0, \
            yi=yinout)

    return (x, yinout)

# y is THE array, holding values that will be integrated in time

N = 1024
dt= 0.00001

x = linspace(0,1,N)
dx= x[1] - x[0]

rhov  = ones(N)
rhov[where((x < 0.7) & (x > 0.4))] = 2.0
rho = ones(N) #ones(N) + peak(N/2, N/6, N, 3)
p_g  = zeros(N) #eos(rho, 5./3)
Q    = zeros(N)
g_z  = 0.0

y    = zeros((4,N))
y[0] = rho
y[1] = rhov
y[2] = p_g
y[3] = Q


fig = plt.figure()
minval = min(x)
maxval = max(x)

# set up axes and plot empty data set
ax  = plt.axes(xlim=(minval, maxval), ylim=(-0.5,3))
line_rho,   = ax.plot([], [], lw=2)
line_rv,    = ax.plot([], [], lw=2)
line_p_g,   = ax.plot([], [], lw=2)
lines = [line_rho, line_rv, line_p_g]

def init():
    for line in lines:
        line.set_data([], [])
    ax.legend([line_rho, line_rv, line_p_g], \
            [r'$\rho$', r'$\rho v$', r'$p_{\rm gas}$'])
    return lines

def animate(it, y, lines):
    x, y = get_state(y)
    rho, rhov, p_g, Q   = y[0], y[1], y[2], y[3]
    lines[0].set_data(x, rho)
    lines[1].set_data(x, rhov)
    lines[2].set_data(x, p_g)
    return lines


anim = animation.FuncAnimation(fig, animate, init_func=init, \
        fargs=(y, lines), frames=64, interval=25, blit=True)

# To save, uncomment following line
# Requires packages ffmpeg, mencoder and libx264 with dependencies
print 'Writing animation to file...'
#anim.save('ec_test.mp4', fps=25, extra_args=['-vcodec', 'libx264', '-profile:v', 'baseline'])
print 'Done'

# the show goes on forever:
plt.show()
