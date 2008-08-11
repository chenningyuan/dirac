#!/usr/bin/env python

import os
import Ngl
import numpy
from Numeric import *

#-------------------------------------------------------------------------------

def set_xyuv_pointer(user_orientation, transpose):

    if user_orientation == 'xy':
       x_pointer = 0
       y_pointer = 1
       u_pointer = 5
       v_pointer = 4

    if user_orientation == 'xz':
       x_pointer = 0
       y_pointer = 2
       u_pointer = 6
       v_pointer = 4

    if user_orientation == 'yz':
       x_pointer = 1
       y_pointer = 2
       u_pointer = 6
       v_pointer = 5

    if transpose:
       save      = x_pointer
       x_pointer = y_pointer
       y_pointer = save
       save      = u_pointer
       u_pointer = v_pointer
       v_pointer = save

    return (x_pointer, y_pointer, u_pointer, v_pointer)

#-------------------------------------------------------------------------------

def plot_figure(figure):

    (x_pointer, y_pointer, u_pointer, v_pointer) = set_xyuv_pointer(figure.user_orientation, figure.transpose)

    dont_want_orientation_l = ['xy', 'xz', 'yz']
    dont_want_orientation_l.remove(figure.user_orientation)

    min_x =  1.0e100
    max_x = -1.0e100
    min_y =  min_x
    max_y =  max_x

    stream_file_s = ''
    scalar_file_s = ''

    if figure.plot_stream:
       f = open(figure.stream_file_name, 'r')
       stream_file_s = f.readlines()
       f.close()

    if figure.plot_scalar:
       f = open(figure.scalar_file_name, 'r')
       scalar_file_s = f.readlines()
       f.close()

    #---------------------------------------------------------------------------

    if figure.plot_stream:

       u_l = []
       v_l = []

       on        = 0
       nr_points = 0

       for line in stream_file_s:
           for orientation in dont_want_orientation_l:
               if orientation in line:
                  on = 0
           if figure.user_orientation in line:
              on = 1
           if on:
              if len(line.split()) == 7:
                 nr_points = nr_points + 1

                 x = float(line.split()[x_pointer])
                 y = float(line.split()[y_pointer])

                 if x < min_x:
                    min_x = x
                 if y < min_y:
                    min_y = y
                 if x > max_x:
                    max_x = x
                 if y > max_y:
                    max_y = y

                 if figure.transpose == 1:
                   u =  float(line.split()[u_pointer])
                   v =  float(line.split()[v_pointer])
                 else:
                   u = -float(line.split()[u_pointer])
                   v = -float(line.split()[v_pointer])

                 u_l.append(u)
                 v_l.append(v)

    #---------------------------------------------------------------------------

    if figure.plot_scalar:

       s_l = []

       on        = 0
       nr_points = 0

       for line in scalar_file_s:
           for orientation in dont_want_orientation_l:
               if orientation in line:
                  on = 0
           if figure.user_orientation in line:
              on = 1
           if on:
              if len(line.split()) == 7:
                 nr_points = nr_points + 1

                 x = float(line.split()[x_pointer])
                 y = float(line.split()[y_pointer])

                 if x < min_x:
                    min_x = x
                 if y < min_y:
                    min_y = y
                 if x > max_x:
                    max_x = x
                 if y > max_y:
                    max_y = y

                 value = figure.scalar_value(float(line.split()[0]), \
                                             float(line.split()[1]), \
                                             float(line.split()[2]), \
                                             float(line.split()[3]), \
                                             float(line.split()[4]), \
                                             float(line.split()[5]), \
                                             float(line.split()[6]))
                 s_l.append(value)

    #---------------------------------------------------------------------------

    nr_points = int(nr_points**0.5)

    if figure.plot_stream:

       u_array = []
       v_array = []

       ipoint = 0
       for i in range(nr_points):

           u_slice = []
           v_slice = []

           for j in range(nr_points):

               norm = (u_l[ipoint]**2.0 + v_l[ipoint]**2.0)**0.5

               u_slice.append(u_l[ipoint]/norm)
               v_slice.append(v_l[ipoint]/norm)

               ipoint = ipoint + 1

           u_array.append(u_slice)
           v_array.append(v_slice)

       u_array = array(u_array, Float0)
       v_array = array(v_array, Float0)

       if figure.transpose == 1:
          u_array = numpy.transpose(u_array)
          v_array = numpy.transpose(v_array)

    #---------------------------------------------------------------------------

    if figure.plot_scalar:

       s_array = []

       ipoint = 0
       for i in range(nr_points):
           s_slice = []

           for j in range(nr_points):
               s_slice.append(s_l[ipoint])
               ipoint = ipoint + 1

           s_array.append(s_slice)

       s_array = array(s_array, Float0)

       if figure.transpose == 1:
          s_array = numpy.transpose(s_array)

    #---------------------------------------------------------------------------

    wks_type = 'ps'
    wks = Ngl.open_wks(wks_type, figure.output_file_name)

    color_res = Ngl.Resources()

    stream_res              = Ngl.Resources()
    stream_res.nglDraw      = False
    stream_res.nglFrame     = False
    stream_res.nglMaximize  = False

    contour_res             = Ngl.Resources()
    contour_res.nglDraw     = False
    contour_res.nglFrame    = False
    contour_res.nglMaximize = False

    #---------------------------------------------------------------------------

    l = []

    l.append(figure.rgb_middle)
    l.append([0.00, 0.00, 0.00])

    f_l = []
    n   = 10

    for f in range(n):
        f_l.append(float((f+1)*(1.0/n)))

    for f in f_l:
        r = figure.rgb_min[0] + f*(figure.rgb_middle[0] - figure.rgb_min[0])
        g = figure.rgb_min[1] + f*(figure.rgb_middle[1] - figure.rgb_min[1])
        b = figure.rgb_min[2] + f*(figure.rgb_middle[2] - figure.rgb_min[2])
        l.append([r, g, b])

    for f in f_l:
        r = figure.rgb_middle[0] + f*(figure.rgb_max[0] - figure.rgb_middle[0])
        g = figure.rgb_middle[1] + f*(figure.rgb_max[1] - figure.rgb_middle[1])
        b = figure.rgb_middle[2] + f*(figure.rgb_max[2] - figure.rgb_middle[2])
        l.append([r, g, b])

    rgb_map_scalar = numpy.array(l, 'f')

    #---------------------------------------------------------------------------

    rgb_min  = [0.9, 0.9, 0.9]
    rgb_max  = [0.0, 0.0, 0.0]

    l = []

    l.append([1.0, 1.0, 1.0])
    l.append([0.0, 0.0, 0.0])

    f_l = []
    n   = 100

    for f in range(n):
        f_l.append(float((f+1)*(1.0/n)))

    for f in f_l:
        r = rgb_min[0] + f*(rgb_max[0] - rgb_min[0])
        g = rgb_min[1] + f*(rgb_max[1] - rgb_min[1])
        b = rgb_min[2] + f*(rgb_max[2] - rgb_min[2])
        l.append([r, g, b])

    rgb_map_stream = numpy.array(l, 'f')

    #---------------------------------------------------------------------------

    stream_res.vpWidthF                 = 0.5
    stream_res.vpHeightF                = 0.5

    stream_res.vcGlyphStyle             = "CurlyVector"
    stream_res.vcMonoFillArrowFillColor = False
    stream_res.vcRefLengthF             = 0.2
    stream_res.vcMinDistanceF           = figure.stream_arrow_min_distance
    stream_res.vcLineArrowThicknessF    = figure.stream_arrow_line_thickness
    stream_res.vcRefAnnoOn              = False

    stream_res.vfXCStartV               = min_x
    stream_res.vfXCEndV                 = max_x
    stream_res.vfYCStartV               = min_y
    stream_res.vfYCEndV                 = max_y

    line_thickness                      = 3.5
    stream_res.tmXBFormat               = "f"
    stream_res.tmYLFormat               = "f"
    stream_res.tmXBMajorThicknessF      = line_thickness
    stream_res.tmXBMinorThicknessF      = line_thickness/2.0
    stream_res.tmXTMajorThicknessF      = line_thickness
    stream_res.tmXTMinorThicknessF      = line_thickness/2.0
    stream_res.tmYLMajorThicknessF      = line_thickness
    stream_res.tmYLMinorThicknessF      = line_thickness/2.0
    stream_res.tmYRMajorThicknessF      = line_thickness
    stream_res.tmYRMinorThicknessF      = line_thickness/2.0
    stream_res.tmBorderThicknessF       = line_thickness

    #---------------------------------------------------------------------------

    contour_res.vpWidthF              = 0.5
    contour_res.vpHeightF             = 0.5

    contour_res.pmLabelBarDisplayMode = "Never"
    contour_res.cnLevelSelectionMode  = "ManualLevels"
    contour_res.cnMinLevelValF        = figure.scalar_range[0]
    contour_res.cnMaxLevelValF        = figure.scalar_range[1]
    contour_res.cnLevelSpacingF       = (figure.scalar_range[1] - figure.scalar_range[0])/20.0

    contour_res.cnLineThicknessF      = figure.contour_line_thickness
    contour_res.cnFillOn              = True
    contour_res.cnMonoFillPattern     = True
    contour_res.cnMonoFillColor       = False
    contour_res.cnLinesOn             = True
    contour_res.cnLineLabelsOn        = False
    contour_res.cnInfoLabelOn         = False

    contour_res.sfXCStartV            = min_x
    contour_res.sfXCEndV              = max_x
    contour_res.sfYCStartV            = min_y
    contour_res.sfYCEndV              = max_y

    contour_res.tmXBFormat            = "f"
    contour_res.tmYLFormat            = "f"

    #---------------------------------------------------------------------------

    if figure.plot_scalar:
       color_res.wkColorMap = rgb_map_scalar
       Ngl.set_values(wks, color_res)
       contour = Ngl.contour(wks, s_array, contour_res)
       Ngl.draw(contour)

    if figure.plot_stream:
       color_res.wkColorMap = rgb_map_stream
       Ngl.set_values(wks, color_res)
       stream = Ngl.vector(wks, u_array, v_array, stream_res)
       Ngl.draw(stream)

    Ngl.frame(wks)
#   Ngl.end()
