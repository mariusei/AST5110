
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
