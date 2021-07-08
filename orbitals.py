import numpy
import scipy.special
import scipy.misc
from enthought.mayavi import mlab

r = lambda x,y,z: numpy.sqrt(x**2+y**2+z**2)
theta = lambda x,y,z: numpy.arccos(z/r(x,y,z))
phi = lambda x,y,z: numpy.arctan(y/x)
#phi = lambda x,y,z: numpy.pi+numpy.select(
#	[x>0, x==0, x<0],
#	[
#		numpy.arctan(y/x),
#		.5*numpy.pi*numpy.sign(y),
#		numpy.arctan(y/x)+numpy.pi*numpy.sign(y)]
#)
a0 = 1.
R = lambda r,n,l: (2*r/n/a0)**l * numpy.exp(-r/n/a0) * scipy.special.genlaguerre(n-l-1,2*l+1)(2*r/n/a0)
WF = lambda r,theta,phi,n,l,m: R(r,n,l) * scipy.special.sph_harm(m,l,phi,theta)
absWF = lambda r,theta,phi,n,l,m: abs(WF(r,theta,phi,n,l,m))**2

x,y,z = numpy.ogrid[-24:24:55j,-24:24:55j,-24:24:55j]

mlab.figure()

#mask = numpy.select([theta(x,y,z)>numpy.pi/3.],[numpy.select([abs(phi(x,y,z))<numpy.pi/3.],[numpy.nan],default=1)],default=1)
mask = 1

for n in range(2,3):
    for l in range(1,n):
        for m in range(-l,l+1,1):
            w = absWF(r(x,y,z),theta(x,y,z),phi(x,y,z),n,l,m)
            mlab.contour3d(w*mask,contours=6,transparent=True)

mlab.colorbar()
mlab.outline()
mlab.show()
