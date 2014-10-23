from stagger import *

N = 1024
x = linspace(0,2*pi,N+1)
dx =x[1]-x[0]

y = sin(x)
x = x[0:N]
y = y[0:N]

yy=deriv(interp(y, direction=-1), dx)
yyy=deriv(interp(yy, direction=-1), dx)
yyyy=deriv(interp(yyy, direction=-1), dx)

y = zeros((4, N))
y[0] = sin(x)
y[1] = yy
y[2] = yyy
y[3] = yyyy

minval = min(x)
maxval = max(x)

#plt.plot(x,y, x, yy, x, yyy, x, yyyy)
#plt.show()


# set up axes and plot empty data set
fig = plt.figure()
ax  = plt.axes(xlim=(minval, maxval), ylim=(-1.5,2.5))
ax.grid()
line_rho,   = ax.plot([], [], lw=2)
line_rv,    = ax.plot([], [], lw=2)
line_p_g,   = ax.plot([], [], lw=2)
line_diffop,   = ax.plot([], [], lw=2)
lines = [line_rho, line_rv, line_p_g, line_diffop]

def get_state(y):
    for i in range(len(y)):
        y[i] = interp(y[i], direction=1)
    return x, y

def init():
    for line in lines:
        line.set_data([], [])
    ax.legend(lines, \
            [r'$y_1=\sin(x)$', r"$y_1'$", r"$y_1''$", r"$y_1'''$"])
    return lines

def animate(it, y, lines):
    x, y = get_state(y)
    rho, rhov, p_g, Q   = y[0], y[1], y[2], y[3]
    lines[0].set_data(x, rho)
    lines[1].set_data(x, rhov)
    lines[2].set_data(x, p_g)
    lines[3].set_data(x, Q)
    return lines


anim = animation.FuncAnimation(fig, animate, init_func=init, \
        fargs=(y, lines), frames=1024, interval=05, blit=True)

# To save, uncomment following line
# Requires packages ffmpeg, mencoder and libx264 with dependencies
print 'Writing animation to file...'
anim.save('stability.mp4', fps=125, extra_args=['-vcodec', 'libx264'])
print 'Done'

# the show goes on forever:
plt.show()
