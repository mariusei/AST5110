from stagger import *

from solvers import *
# Solvers take functions: 
#   f(x_in, y_in, dx, dt)
# as argument F in
#   solver(F, dx, dt, x_in, y_in)


def equation_set(xi, yi, dx, dt):

    pin  = yi[0]
    rhoi = yi[1]

    Q    = zeros(N)
    g_z  = 0.0
    gamma= 5./3

    dp_dt =     - deriv( 1./rhoi * interp( pin**2, direction=-1 ),  dx) \
                + deriv( eos(rhoi, gamma) + Q,   dx) \
                + interp( rhoi, direction=+1 ) * g_z

    drho_dt =   - deriv( pin, dx, direction=-1 )

    return array([dp_dt, drho_dt])

def eos(rho, gamma):
    rho[where(rho <= 0)] = 0
    return power(rho, gamma)


## INITIAL CONDITIONS ##

N = 1024
x = linspace(0,10,N)
y = zeros((2,N))

dx = x[0] - x[1]
dt = dx/5

# Momentum
y[0,:] = 0

# Density
y[1, :]                         = 0.8
y[1, where((x < 6) & (x > 3))]  = 1.0



## ANIMATION ##

def get_state(yinout):

    yinout[:] = rk4solver(f=equation_set, dx=dx, dt=dt, xi=zeros(N), \
            yi=yinout)

    return (x, yinout)

fig = plt.figure()
minval = min(x)
maxval = max(x)

# set up axes and plot empty data set
ax  = plt.axes(xlim=(minval, maxval), ylim=(-0.8,01.5))
line_rhov,  = ax.plot([], [], lw=2)
line_rho,   = ax.plot([], [], lw=2)
lines = [line_rhov, line_rho]

def init():
    for line in lines:
        line.set_data([], [])
    ax.legend(lines, \
            [r'$\rho v$', r'$\rho$'])
    return lines

def animate(it, y, lines):
    x, y = get_state(y)
    rhov  = y[0]
    rho   = y[1]
    lines[0].set_data(x, rhov)
    lines[1].set_data(x, rho)
    return lines


anim = animation.FuncAnimation(fig, animate, init_func=init, \
        fargs=(y, lines), frames=512, interval=205, blit=True)

# To save, uncomment following line
# Requires packages ffmpeg, mencoder and libx264 with dependencies
print 'Writing animation to file...'
#anim.save('working.mp4', fps=240, extra_args=['-vcodec', 'libx264'])
print 'Done'

# the show goes on forever:
plt.show()




