#!/usr/bin/env python
#
# script to convert raw 2D plots generated with DIRAC
# to nice and shiny ps files using PyNGL
#
# requires: http://www.pyngl.ucar.edu/
#           http://numpy.scipy.org/
#
# you can plot a streamline plot or a scalar contour plot
# or both in one figure
#
# DIRAC 2D format:
# each line contains seven numbers: x, y, z, s, u, v, w
#
#                             with: x, y, z  position
#                                   s        scalar value
#                                   u, v, w  vector components
#
# with these seven numbers you can form other scalar functions
# using the function scalar_value(x, y, z, s, u, v, w)
#
# (radovan bast 2008-08-11)

from visual_common import *

class figure_c:

      def __init__(self):

          self.transpose = 1

          self.plot_stream = 1
          self.plot_scalar = 1

          self.user_orientation = 'xy'          # 'xy' or 'xz' or 'yz'

          self.scalar_range = [-1.0e+1, 1.0e+1] # this will be used to adjust coloring of the scalar plot

          self.rgb_min    = [0.0, 0.0, 1.0]     # color for scalar_range[0]
          self.rgb_middle = [1.0, 1.0, 1.0]
          self.rgb_max    = [1.0, 0.0, 0.0]     # color for scalar_range[1]
                                                # so color it from blue (minus) over white to red (plus)
          self.stream_file_name = ''
          self.scalar_file_name = ''
          self.output_file_name = ''

          self.stream_arrow_line_thickness = 0.5
          self.stream_arrow_min_distance   = 0.01
          self.contour_line_thickness      = 0.1

      def scalar_value(self, x, y, z, s, u, v, w):
          scalar_value = (u**2.0 + v**2.0 + w**2.0)**0.5
          return scalar_value

#-------------------------------------------------------------------------------

example_figure                  = figure_c()
example_figure.stream_file_name = 'example.2d'
example_figure.plot_scalar      = 0
example_figure.output_file_name = 'example1'

plot_figure(example_figure)


example_figure                  = figure_c()
example_figure.plot_stream      = 0
example_figure.scalar_file_name = 'example.2d'
example_figure.output_file_name = 'example2'

plot_figure(example_figure)


example_figure                  = figure_c()
example_figure.stream_file_name = 'example.2d'
example_figure.scalar_file_name = 'example.2d'
example_figure.output_file_name = 'example3'

plot_figure(example_figure)
