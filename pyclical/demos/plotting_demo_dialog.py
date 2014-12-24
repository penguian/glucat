#!/usr/bin/env python
"""
An example of how to modify the data visualized  via an interactive dialog.
"""
# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# Copyright (c) 2008, Enthought, Inc.
# License: BSD Style.

from traits.api import Bool, Button, HasTraits, Instance, Range, on_trait_change
from traitsui.api import Group, View

from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MlabSceneModel

from plotting_demo_mayavi import demo,  \
        default_nbr_orbits,             \
        default_nbr_points,             \
        default_scaling,                \
        default_segment_len,            \
        default_figwidth,               \
        default_figheight,              \
        default_radius,                 \
        default_resolution,             \
        default_opacity

class TorusDemoDialog(HasTraits):
    """
    Parameters:
    number_of_orbits    : Number of orbits to plot.
    number_of_points    : Number of points to plot in each orbit.
    scaling             : Scaling constant to use with bivector br.
    segment_length      : Number of points in an orbit segment.
    figure_width        : Width of figure.
    figure_height       : Height of figure.
    radius              : Relative radius of the ball around each point.
    resolution          : Resolution of the ball around each point.
    opacity             : Opacity of the ball around each point.
    reciprocal          : Use the reciprocal bivector with respect to the pseudoscalar of R_{4,0}
    """
    number_of_orbits    = Range(  1,     16, default_nbr_orbits)
    number_of_points    = Range(  2,  32000, default_nbr_points)
    scaling             = Range(  1,  32000, default_scaling)
    segment_length      = Range(  1,  32000, default_segment_len)
    figure_width        = Range(128,   4096, default_figwidth)
    figure_height       = Range(128,   4096, default_figwidth)
    radius              = Range(  0.01,   0.1, default_radius)
    resolution          = Range(  3,     32, default_resolution)
    opacity             = Range(  0.0,    1.0, default_opacity)
    use_reciprocal      = Bool
    plot_button         = Button('Draw a new plot')

    scene = Instance(MlabSceneModel, ())

    plot = Instance(PipelineBase)


    # When the scene is activated, or when the parameters are changed, we
    # update the plot.
    @on_trait_change('plot_button')
    def new_plot(self):
        demo(
            nbr_figures = 1,
            nbr_orbits  = self.number_of_orbits,
            nbr_points  = self.number_of_points,
            scaling     = self.scaling,
            segment_len = self.segment_length,
            figwidth    = self.figure_width,
            figheight   = self.figure_height,
            radius      = self.radius,
            resolution  = self.resolution,
            opacity     = self.opacity,
            reciprocal  = self.use_reciprocal)

    # The layout of the dialog created
    view = View(Group(
                    '_',
                    'number_of_orbits',
                    'resolution',
                    'use_reciprocal',
                    '_',
                    'plot_button',
                    label='Plotting controls',
                    ),
                Group(
                    '_',
                    'number_of_points',
                    'segment_length',
                    'scaling',
                    'figure_width',
                    'figure_height',
                    label='Advanced controls',
                    ),
                resizable=True,
                title='Torus demo'
                )

this_dialog = TorusDemoDialog()
this_dialog.configure_traits()
