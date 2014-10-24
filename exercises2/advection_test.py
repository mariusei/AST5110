
from stagger import *

def ecsolver(f, dx, dt, xi, yi):
    """ Advances one step in time """
    y1 = yi + dt * f(xi, yi, dx, dt)

    return y1

def rk4solver(f, dx, dt, xi, yi):
    """ Advances one step in time using a predictor-corrector method """

    k1 = f(xi, yi, dx, dt)
    k2 = f(xi, yi + dt/2 * k1, dx, dt/2)
    k3 = f(xi, yi + dt/2 * k2, dx, dt/2)
    k4 = f(xi, yi + dt   * k3, dx, dt)

    y1 = yi + dt/6 * (k1 + 2.*k2 + 2.*k3 + k4)

    return y1

def advection_eqs(xi, yi, dx, dt):
    """ Advection equation from first exercise set

    d rho/dt = -v d rho/dz

    """

    v   = yi[0]
    rho = yi[1]
    rho[where(rho < 0)] = 0.001
    
    # To find diffusion operator
    gamma = 5./3
    K     = 1.0

    p_g  = K * eos(rho, gamma)

    nu   = 1.0  # of order unity
    c_s  = K * gamma * power( rho, gamma - 1)
    mu   = rho * c_s * nu * dx

    Q    =      mu * deriv(v, dx)



    dvdt    = zeros(N)

    drhodt  = - interp(v, direction=+1) * deriv( \
                                interp(rho,direction=1), \
                                dx, direction=-1) \
              + Q

    return array([dvdt, drhodt, Q])

def eos(rho, gamma):
    return power(rho, gamma)

# Initialise

N = 1024
x = linspace(0,1,N)
dx = x[1] - x[0]
dt = dx/10

y = zeros((3,N))

y[0,:] = -15.7
y[1,:] = 0.1
y[1, (where((x > 0.3) & (x < 0.4)))] = 1.0


### ANIMATION ###


def get_state(yinout):

    yinout[:] = rk4solver(f=advection_eqs, dx=dx, dt=dt, xi=0, \
            yi=yinout)

    return (x, yinout)

fig = plt.figure()
minval = min(x)
maxval = max(x)

# set up axes and plot empty data set
ax  = plt.axes(xlim=(minval, maxval), ylim=(-2.5,03))
line_rho,   = ax.plot([], [], lw=2)
line_Q,   = ax.plot([], [], lw=2)
lines = [line_rho, line_Q]

def init():
    for line in lines:
        line.set_data([], [])
    ax.legend(lines, \
            [r'$\rho$', r'Q'])
    return lines

def animate(it, y, lines):
    x, y = get_state(y)
    rho   = y[1]
    Q     = y[2]
    lines[0].set_data(x, rho)
    lines[1].set_data(x, Q)
    return lines


anim = animation.FuncAnimation(fig, animate, init_func=init, \
        fargs=(y, lines), frames=512, interval=005, blit=True)

# To save, uncomment following line
# Requires packages ffmpeg, mencoder and libx264 with dependencies
print 'Writing animation to file...'
#anim.save('damped.mp4', fps=50, extra_args=['-vcodec', 'libx264'])
print 'Done'

# the show goes on forever:
plt.show()
