import numpy as np
import matplotlib.pyplot as plt

def rk4(f, x, t, dt, stages=4, s=0):
    """Runge-Kutta (explicit, non-adaptive) numerical (S)ODE solvers.

    For ODEs, the order of convergence equals the number of `stages`.

    For SDEs with additive noise (`s>0`), the order of convergence
    (both weak and strong) is 1 for `stages` equal to 1 or 4.
    These correspond to the classic Euler-Maruyama scheme and the Runge-Kutta
    scheme for S-ODEs respectively, see `bib.grudzien2020numerical`
    for a DA-specific discussion on integration schemes and their discretization errors.

    Parameters
    ----------
    f : function
        The time derivative of the dynamical system. Must be of the form `f(t, x)`

    x : ndarray or float
        State vector of the forcing term

    t : float
        Starting time of the integration

    dt : float
        Integration time step.

    stages : int, optional
        The number of stages of the RK method.
        When `stages=1`, this becomes the Euler (-Maruyama) scheme.
        Default: 4.

    s : float
        The diffusion coeffient (std. dev) for models with additive noise.
        Default: 0, yielding deterministic integration.

    Returns
    -------
    ndarray
        State vector at the new time, `t+dt`
    """

    # Draw noise
    if s > 0:
        W = s * np.sqrt(dt) * np.random.randn(*x.shape)
    else:
        W = 0

    # Approximations to Delta x
    if stages >= 1: k1 = dt * f(x,           t)         + W    # noqa
    if stages >= 2: k2 = dt * f(x+k1/2.0,    t+dt/2.0)  + W    # noqa
    if stages == 3: k3 = dt * f(x+k2*2.0-k1, t+dt)      + W    # noqa
    if stages == 4:
                    k3 = dt * f(x+k2/2.0,    t+dt/2.0)  + W    # noqa
                    k4 = dt * f(x+k3,        t+dt)      + W    # noqa

    # Mix proxies
    if    stages == 1: y = x + k1                              # noqa
    elif  stages == 2: y = x + k2                              # noqa
    elif  stages == 3: y = x + (k1 + 4.0*k2 + k3)/6.0          # noqa
    elif  stages == 4: y = x + (k1 + 2.0*(k2 + k3) + k4)/6.0   # noqa
    else:
        raise NotImplementedError

    return y


def x0(M):
    x = np.zeros(M)
    x[0] = 1
    return x


def shift(x, n):
    return np.roll(x, -n, axis=-1)


def dxdt_autonomous(x):
    return (shift(x, 1)-shift(x, -2))*shift(x, -1) - x


def dxdt(x):
    return dxdt_autonomous(x) + 8


def step(x0, t, dt):
    return rk4(lambda x, t: dxdt(x), x0, np.nan, dt)

def plot():
    theta = np.linspace(0, 2*np.pi, nx)
    fig = plt.figure(1)
    ax = fig.add_subplot(polar=True)
    for it in range(nt):
        plt.cla()
        ax.plot(theta, x[it])
        ax.set_rmin(-10)
        ax.set_rmax(10)
        ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
        ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
        ax.grid(True)

        plt.pause(0.05)
    plt.show()

if __name__ == "__main__":

    nt = 20000
    nx = 6
    dt = 0.05
    # initial condition
    x_spinup = np.random.random(nx)
    # truth trajectory
    x = np.zeros((nt + 1, nx))
    # obs trajectory
    y = np.zeros((nt + 1, nx))
    for it in range(nt):
        x_spinup = step(x_spinup, dt*it, dt)

    x[0] = x_spinup
    y[0] = x[0] + np.random.randn(nx)

    for it in range(nt):
        x[it + 1] = step(x[it], dt*it, dt)
        y[it + 1] = x[it + 1] + np.random.randn(nx)

    np.savez_compressed('trajectory.npz', truth = x, obs = y)


