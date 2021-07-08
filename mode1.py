import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# The following import configures Matplotlib for 3D plotting.
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm, genlaguerre
plt.rc('text', usetex=True)

# Grids of polar and azimuthal angles
theta = np.linspace(0, np.pi, 300)
phi = np.linspace(0, 2*np.pi, 300)
# Create a 2-D meshgrid of (theta, phi) angles.
theta, phi = np.meshgrid(theta, phi)
# Calculate the Cartesian coordinates of each point in the mesh.
xyz = np.array([np.sin(theta) * np.sin(phi),
                np.sin(theta) * np.cos(phi),
                np.cos(theta)])

a0 = 1.0
r0 = 2.0
R = lambda r, n, l: (2*r/n/a0)**l * np.exp(-r/n/a0) * genlaguerre(n-l-1,2*l+1)(2*r/n/a0)

def plot_Y(ax, el, m, n):
    """Plot the spherical harmonic of degree el and order m on Axes ax."""

    # NB In SciPy's sph_harm function the azimuthal coordinate, theta,
    # comes before the polar coordinate, phi.
    Y = sph_harm(abs(m), el, phi, theta)

    # Linear combination of Y_l,m and Y_l,-m to create the real form.
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    # Yx, Yy, Yz = np.abs(Y) * xyz

    Rs = np.linalg.norm(Y, axis=1)

    WF = R(Rs, n, l) * Y
    # Yx, Yy, Yz = np.abs(WF) * xyz
    Yx, Yy, Yz = np.abs(WF)**2 * xyz


    # Colour the plotted surface according to the sign of Y.
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('binary'))
    cmap.set_clim(0, .95)

    ax.plot_surface(Yx, Yy, Yz,
                    facecolors=cmap.to_rgba(WF.real),
                    rstride=2, cstride=2,
                    antialiased=True)

    # Set the Axes limits and title, turn off the Axes frame.
    # ax.set_title(r'$Y_{{{},{}}}$'.format(el, m))
    ax_lim = .5
    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.axis('on')

fig = plt.figure(figsize=(8, 8), dpi=100)
fig.tight_layout(pad=0.0)

def plotblablabla(l, m, n):
    print('n: {}, l: {}, m: {}'.format(n, l, m))
    ax = fig.add_subplot(projection='3d')
    plot_Y(ax, l, m, n)
    plt.savefig('Y{}_{}_{}.png'.format(n, l, m), transparent=True, bbox_inches='tight')
    fig.clear()

for n in range(1, 3):
    for l in range(0, n):
        for m in range(-l, l + 1):
            plotblablabla(l, m, n)
