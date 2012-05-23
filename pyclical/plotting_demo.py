# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# plotting_demo.py: Demonstrate the use of Matplotlib plotting with PyClical.
#
#    copyright            : (C) 2010-2012 by Paul C. Leopardi
#
#    This library is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library.  If not, see <http://www.gnu.org/licenses/>.

from PyClical import *

#
# Default values common to functions draw_orbit and demo.
#
default_nbr_points  = 40000
default_segment_len =  8000
default_azimuth     =   210.0
default_rot_angle   =   360
default_jitter      =    30

def draw_orbit(r,
        nbr_points  = default_nbr_points,
        segment_len = default_segment_len,
        azimuth     = default_azimuth,
        rot_angle   = default_rot_angle,
        jitter      = default_jitter):
    """
    Plot a curve created by the rotor r,

    Parameters:
    r           : Rotor, an even element of R_{4,1}.
                  assumed to be an element of Spin(4,1).
    nbr_points  : Number of points overall.
    segment_len : Number of points in a curve segment.
    azimuth     : Angle about a vertical axis used to define the viewpoint.
    rot_angle   : Angle in degrees to use to rotate each curve for display.
    jitter      : Angle in degrees to use for each step in the rotation.
    """
    #
    # Imports needed for array calculation and plotting.
    #
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    #
    # Array of points.
    #
    p = np.empty((segment_len, 3))
    #
    # Use a new figure.
    #
    fig = plt.figure(figsize = (15, 12))
    ax  = fig.gca(projection = '3d')
    ax.view_init(azim=azimuth)
    #
    # Draw the origin as a white point. This works around a bug in mplot3d where
    # plots with a small range in Z values are plotted as if Z is close to 0.
    #
    ax.scatter([0], [0], [0], c='white', alpha=0.5, edgecolors='none')
    plt.draw()
    #
    # Initial CGA null vector corresponding to point
    # on the unit sphere in 3D Euclidean space.
    #
    x = random_clifford(istpq(3,0))(1)
    x = e(1) if x == 0 else x / abs(x)
    u = cga3(x)(1)
    #
    # Find the inverse of the rotor r.
    #
    ir = inv(r)
    #
    # Number of segments.
    #
    M = nbr_points / segment_len
    #
    # Split the curve into M segments.
    #
    for j in xrange(M):
        #
        # Find segment_len points forming a curve segment
        # by successively using the rotor r and its inverse.
        #
        for k in xrange(segment_len):
            #
            # Determine the current 3D Euclidean point
            # corresponding to the CGA null vector u.
            #
            p[k, :] = agc3(u).vector_part()
            #
            # Act on u via the adjoint action of the rotor r.
            #
            u = r * u * ir
        #
        # Calculate the norms of the points in p and store them in n.
        #
        n = np.sqrt(np.reshape(np.sum(p*p, 1), [segment_len, 1]))
        #
        # Fudge norm 0 to the value 1 to avoid dividing by 0.
        #
        n[n==0] = 1.0
        #
        # Tile the norm n to be a segmen_len by 3 array
        # for use in determining colours.
        #
        n3 = np.tile(n, (1, 3))
        #
        # Set c to be an array of vectors in the unit ball with directions
        # corresponding to p. Use c to determine colours.
        # Try to make colours bright without oversaturating them:
        # arctan(x)/2 goes to pi/4 < 1 as x goes to infinity.
        #
        c = (p / n3) * (np.arctan(n3) / 2.0)
        #
        # Plot the new scattered points using the colours given by (c+1)/2.
        #
        ax.scatter(p[:, 0], p[:, 1], p[:, 2], c=(c + 1.0) / 2.0, alpha=0.5, edgecolors='none')
        plt.draw()
    #
    # Rotate the plot about a vertical axis by rot_angle degrees
    #
    for phi in range(jitter, rot_angle + jitter, jitter):
        ax.view_init(azim=azimuth + phi)
        plt.draw()

#
# Default values for demo.
#
default_nbr_curves  =    10
default_scaling     =  4000

def demo(d=4,
        nbr_curves  = default_nbr_curves,
        nbr_points  = default_nbr_points,
        scaling     = default_scaling,
        segment_len = default_segment_len,
        azimuth     = default_azimuth,
        rot_angle   = default_rot_angle,
        jitter=default_jitter):
    """
    Plot curves created by exponentiating a random bivector in R_{d,0}.

    Parameters:
    d           : Dimension of bivector.
    nbr_curves  : Number of curves to plot.
    nbr_points  : Number of points overall.
    scaling     : Scaling constant to use with bivector.
    segment_len : Number of points in a curve segment.
    azimuth     : Angle about a vertical axis used to define the viewpoint.
    rot_angle   : Angle in degrees to use to rotate each curve for display.
    jitter      : Angle in degrees to use for each step in the rotation.
    """
    #
    # Plot nbr_curves curves.
    #
    for i in xrange(nbr_curves):
        #
        # Use a random bivector in R_{d,0} with appropriate scaling.
        #
        b = random_clifford(istpq(d,0))(2) * scaling / nbr_points
        #
        # Exponentiate the bivector to obtain a rotor.
        #
        r  = exp(b)
        #
        # Draw the curve.
        #
        draw_orbit(r, nbr_points, segment_len, azimuth, rot_angle, jitter)
