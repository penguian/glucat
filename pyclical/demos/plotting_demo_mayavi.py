#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# plotting_demo_mayavi.py: Demonstrate the use of Mayavi2 plotting with PyClical.
#
#    References:
#    [B] Michael F. Barnsley, Superfractals, http://www.superfractals.com/
#    [DV] Leo Dorst and Robert Valkenburg, "Square Root and Logarithm of Rotors in 3D
#    Geometric Algebra Using Polar Decomposition", in Leo Dorst and Joan Lasenby
#    (eds.), Guide to Geometric Algebra in Practice, Springer, 2011, pp. 81-104.
#
#    copyright            : (C) 2010-2014 by Paul C. Leopardi
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

#
# Imports needed for array calculation and plotting.
#
from builtins import range
import numpy as np
import mayavi.mlab as ml
from PyClical import *

#
# Default values common to functions draw_orbit and demo.
#
default_nbr_points  =  8000
default_segment_len =  8000
default_figwidth    =  1024
default_figheight   =  1024
default_radius      =     0.05
default_resolution  =     8
default_opacity     =     1.0

def draw_orbit(r, s,
        nbr_points  = default_nbr_points,
        segment_len = default_segment_len,
        figwidth    = default_figwidth,
        figheight   = default_figheight,
        radius      = default_radius,
        resolution  = default_resolution,
        opacity     = default_opacity):
    """
    Plot an orbit created by a random sequence using the rotors r and s,

    Parameters:
    r, s        : Rotors, even elements of R_{4,1}.
                  assumed to be elements of Spin(4,1).
    nbr_points  : Number of points overall.
    segment_len : Number of points in an orbit segment.
    figwidth    : Width of figure.
    figheight   : Height of figure.
    radius      : Relative radius of the ball around each point.
    resolution  : Resolution of the ball around each point.
    opacity     : Opacity of the ball around each point.
    """
    #
    # Frame for 3D Euclidean space R^3.
    #
    r3frame = istpq(3,0)
    #
    # Frame for Conformal Geometric Algebra (CGA).
    #
    cga3frame = istpq(4,1)
    #
    # Reframe the rotors r and s, for speed.
    #
    r = r.reframe(cga3frame)
    s = s.reframe(cga3frame)
    #
    # Initial CGA null vector corresponding to point
    # on the unit sphere in 3D Euclidean space.
    #
    x = random_clifford(r3frame)(1)
    x = e(1) if x == 0 else x / abs(x)
    u = cga3(x)(1)
    #
    # Find the inverses of the rotors r and s.
    #
    invr = inv(r)
    invs = inv(s)
    #
    # Array of points to use in plotting.
    #
    p = np.empty((segment_len, 3))
    #
    # Split the orbit into M segments.
    #
    M = nbr_points // segment_len
    for j in range(M):
        #
        # Find segment_len points forming an orbit segment
        # by successively using the rotor r and its inverse.
        #
        for k in range(segment_len):
            #
            # Determine the current 3D Euclidean point
            # corresponding to the CGA null vector u.
            #
            p[k,:] = agc3(u).vector_part(r3frame)
            #
            # Act on u via the adjoint action of either the rotor r or the rotor s.
            # See [DV] for orbits of rotors in 3D CGA.
            # The rotor is chosen uniformly at random, making this a simple chaos game.
            # See [B] for related ideas.
            #
            if np.random.rand() < 0.5:
                u = r * u * invr
            else:
                u = s * u * invs
        #
        # Calculate the norms of the points in p and store them in n.
        #
        n = np.sqrt(np.reshape(np.sum(p*p, 1), [segment_len, 1]))
        #
        # Plot the new scattered points.
        #
        ml.points3d(p[:,0], p[:,1], p[:,2], n[:,0],
                    colormap="jet",
                    scale_factor=radius,
                    resolution=resolution,
                    opacity=opacity)
    #
    # Show the plot.
    #
    ml.show()

#
# Default values for demo.
#
default_nbr_figures =    4
default_nbr_orbits  =    1
default_scaling     = 1000
default_reciprocal  = True

def demo(
        nbr_figures = default_nbr_figures,
        nbr_orbits  = default_nbr_orbits,
        nbr_points  = default_nbr_points,
        scaling     = default_scaling,
        segment_len = default_segment_len,
        figwidth    = default_figwidth,
        figheight   = default_figheight,
        radius      = default_radius,
        resolution  = default_resolution,
        opacity     = default_opacity,
        reciprocal  = default_reciprocal):
    """
    Plot orbits created by exponentiating a random bivector and its reciprocal in R_{4,0}.

    Parameters:
    nbr_figures : Number of figures to plot.
    nbr_orbits  : Number of orbits to plot per figure.
    nbr_points  : Number of points to plot in each orbit.
    scaling     : Scaling constant to use with bivector br.
    segment_len : Number of points in an orbit segment.
    figwidth    : Width of figure.
    figheight   : Height of figure.
    radius      : Relative radius of the ball around each point.
    resolution  : Resolution of the ball around each point.
    opacity     : Opacity of the ball around each point.
    reciprocal  : Use the reciprocal bivector with respect to the pseudoscalar of R_{4,0}
    """
    #
    # Check and adjust arguments for sanity.
    #
    if segment_len > nbr_points:
        segment_len = nbr_points
    #
    # Frame for 4D Euclidean space R^4.
    #
    r4frame = istpq(4,0)
    #
    # Plot nbr_figures figures.
    #
    for fig in range(nbr_figures):
        #
        # Use a new figure.
        #
        ml.figure(size=(figwidth, figheight))
        #
        # Plot nbr_orbits orbits.
        #
        for i in range(nbr_orbits):
            #
            # Set br to be a random bivector in R_{4,0} with appropriate scaling.
            #
            br = random_clifford(r4frame)(2) * scaling / nbr_points
            #
            # Set bs to be the reciprocal bivector with respect to the pseudoscalar of R_{4,0}.
            #
            bs = cl(r4frame) / br if reciprocal else cl(0)
            #
            # Exponentiate the bivectors br and bs to obtain rotors r and s.
            #
            r  = exp(br)
            s  = exp(bs)
            #
            # Draw the orbit.
            #
            draw_orbit(r, s,
                       nbr_points=nbr_points,
                       segment_len=segment_len,
                       figwidth=figwidth,
                       figheight=figheight,
                       radius=radius,
                       resolution=resolution,
                       opacity=opacity)


if __name__ == "__main__":
    demo()
