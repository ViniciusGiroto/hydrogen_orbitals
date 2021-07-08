import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# The following import configures Matplotlib for 3D plotting.
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm, genlaguerre, factorial, lpmv
plt.rc('text', usetex=True)

a0 = 0.5292
points = 400

def L(l, m, rho):
    sum = 0
    for k in range(0, n - l):
        arr = [
            (-1.0) ** (k + 1),
            factorial(n + l) ** 2.0,
            factorial(n - l - 1 - k),
            factorial(2 * l + 1 + k),
            factorial(k),
            rho ** k,
        ]
        sum += arr[0]*arr[1]/(arr[2]*arr[3]*arr[4])*arr[5]

    return sum

def R(r, n, l):
    rho = 2*r/(n*a0)
    sum = L(l, m, rho)

    arr = [
        (2.0 / (n * a0)) ** 3.0,
        factorial(n - l - 1),
        2.0 * n,
        factorial(n + l) ** 3.0,
    ]

    normTermSquared = arr[0]*arr[1]/(arr[2]*arr[3])
    radialTerm = np.exp(-rho / 2.0) * (rho ** l)

    return (sum * radialTerm) ** 2 * normTermSquared

def angDF(l, m, theta, phi):
    absm = np.abs(m)
    normRoot = factorial(l - absm)/factorial(l + absm)
    normRoot = normRoot * (2.0 * l + 1.0)/(4.0 * np.pi)
    legendre = lpmv(absm, l, np.cos(theta))

    return normRoot * legendre ** 2.0 * np.exp(m * 2.0j * phi)


def plot_Y(fig, ax, n, l, m, scale):
    """Plot the spherical harmonic of degree el and order m on Axes ax."""
    r0 = scale

    x = np.linspace(-r0, r0, points)
    y = np.linspace(r0, -r0, points)

    WF = np.zeros((points, points))
    s = 0

    for dx in range(points):
        for dy in range(points):
            kx = x[dx]
            ky = y[dy]
            phi = math.atan2(ky, kx)
            theta = math.atan2(np.hypot(kx, ky), 0)
            r = np.hypot(kx, ky)
            # theta e phi ao contrário
            # Y = angDF(l, m, phi, theta)
            Y = sph_harm(m, l, theta, phi)

            # if m < 0:
            #     Y = np.sqrt(2) * (-1)**m * Y.imag
            # elif m > 0:
            #     Y = np.sqrt(2) * (-1)**m * Y.real

            # k = np.abs(R(r, n, l) * Y)**2
            # wf = Y * R(r, n, l) # R(r, n, l) * Y ** 2
            wf = Y ** 2 * R(r, n, l)
            prob = np.abs(wf)
            s += prob

            WF.itemset((dx, dy), prob)
            # WF[dx, dy] = np.abs(R(r, n, l) * Y)**2

    # WF = np.log10(WF + 1)
    # cmap.set_clim(0, .5)

    # plot R
    """
    kx = np.linspace(0, r0, 100)
    ky = list(map(lambda r: R(r, n, l), kx))

    ax.plot(kx, ky)
    ax.set_xticks(np.linspace(0, r0, 6))
    ax.set_xticklabels(range(0, 6))
    """

    # plot bla
    min = np.min(WF)
    max = np.max(WF)
    print('min: {}, max: {}, sum: {}'.format(min, max, s))

    # plt.grid(True)
    im = ax.imshow(WF, cmap='binary')
    ax.tick_params(direction='in')
    ax.set_xticks(np.linspace(0, points, 5))
    ax.set_yticks(np.linspace(0, points, 5))
    ax.set_xticklabels(map(lambda x: r'$\mathbf{{{}}}$'.format(int(x)), np.linspace(-scale, scale, 5)))
    ax.set_yticklabels(map(lambda x: r'$\mathbf{{{}}}$'.format(int(x)), np.linspace(scale, -scale, 5)))
    ax.set_xlabel(r'$\bm{x \left(\mathrm{\AA}\right)}$')
    ax.set_ylabel(r'$\bm{y \left(\mathrm{\AA}\right)}$')
    ax.set_title(r'$\bm{{({}, {}, {})}}$'.format(n, l, m))
    plt.setp(ax.spines.values(), linewidth=2)
    # fig.colorbar(im, ax=ax)

    # Colour the plotted surface according to the sign of Y.

def plotblablabla(n, l, m, scale):
    fig = plt.figure(figsize=(2.5, 2.5), dpi=400)
    plt.rc('text.latex', preamble=r'\usepackage{bm}')
    fig.tight_layout(pad=0.0)
    ax = plt.axes()
    plot_Y(fig, ax, n, l, m, scale)
    plt.savefig('Y_{}_{}_{}.png'.format(n, l, m), transparent=False, bbox_inches='tight')
    plt.close()

bla = False

n, l, m = (10, 8, 2)
plotblablabla(10, 8, 2, 150)

if bla:
    for n in range(1, 7):
        for l in range(0, n):
            for m in range(-l, l + 1):
                print('n: {}, l: {}, m: {}'.format(n, l, m))
                plotblablabla(n, l, m, 20)
