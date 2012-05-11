# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# plotting_demo.py: Demonstrate the use of Matloblib plotting with PyClical.
#
#    copyright            : (C) 2008-2012 by Paul C. Leopardi
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

def demo(d=4, nbr_curves=2, nbr_points=10000, scaling=200, segment_len=100, arena_width=10.0, rot_angle=30):
    """
    Plot curves created by exponentiating a random bivector in R_{d,0}.
    
    Parameters:
    d           : Dimension of bivector.
    nbr_curves  : Number of curves to plot.
    nbr_points  : Number of points overall.
    scaling     : Scaling constant to use with bivector.
    segment_len : Number of points in a curve segment.
    arena_width : Width of region used to determine colours for scatter plot.
    rot_angle   : Angle in degrees to use to rotate each curve for display.
    """
    #
    # Constant to control the plotting.
    #
    M = nbr_points / segment_len
    #
    # Imports needed for array calculation and plotting.
    #
    from numpy import empty, array, mod
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    #
    # Array of points.
    #
    p = empty((3, nbr_points))
    #
    # Array of RGB colour values.
    # 
    rgb = empty((nbr_points, 3))
    #
    # Plotting constants
    #
    azimuth = 30
    half_arena_width = arena_width / 2.0
    #
    # Plot nbr_curves curves.
    #
    for i in xrange(nbr_curves):
        #
        # Use a new figure for each curve.
        #
        fig = plt.figure(figsize = (15, 12))
        ax  =  fig.gca(projection = '3d')
        ax.view_init(azim=azimuth)
        plt.show()
        #
        # Initial point.
        #
        u = cga3(random_clifford(istpq(3,0))(1))
        #
        # Use a random bivector in R_{d,0} with appropriate scaling.
        #
        b = random_clifford(istpq(d,0))(2) * scaling / nbr_points
        #
        # Exponentiate the bivector to obtain a rotor.
        #
        r = exp(b)
        #
        # Split the curve into M segments.
        #
        for j in xrange(M):
            #
            # Find segment_len points forming a curve segment
            # by successively using the rotor r.
            #
            abot = j*segment_len
            atop = abot + segment_len
            for k in xrange(abot, atop):
                p[:, k] = agc3(u).vector_part()
                u |= r
            #
            # Plot the curve segment.
            #
            alow = abot-1 if j > 0 else abot
            X = p[0, alow:atop]
            Y = p[1, alow:atop]
            Z = p[2, alow:atop]
            ax.plot(X, Y, Z, c='gray', linestyle=':')
            #
            # Determine the colour of the scattered points.
            #
            for h in xrange(3):
                rgb[abot:atop, h] = mod(p[h, abot:atop] + half_arena_width, arena_width) / arena_width
            #
            # Plot the scattered points using the chosen colour.
            #
            ax.scatter(X, Y, Z, c=rgb[alow:atop, :], alpha=0.5, edgecolors='none')
            plt.draw()
        #
        # Rotate the plot about a vertical axis by rot_angle degrees
        #
        for phi in xrange(rot_angle):
            ax.view_init(azim=azimuth + phi)
            plt.draw()
