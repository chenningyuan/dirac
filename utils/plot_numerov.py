# requires PyNGL http://www.pyngl.ucar.edu/
#
# last modification: 25/01/2008 (radovan bast)

import Ngl
from pylab import *

def plot_numerov(potential_x_l,
                 potential_y_l,
                 potential_coefficients,
                 property_x_l,
                 property_y_l,
                 property_coefficients,
                 q,
                 psi_squared_save,
                 energy_save,
                 x_name,
                 potential_name,
                 property_name,
                 file_name):

    wks = Ngl.open_wks('ps', file_name)

    res = Ngl.Resources()

    res.nglDraw            = False
    res.nglFrame           = False
    res.nglMaximize        = False
    res.tiMainString       = ''

    res.trXMaxF            = 0.5
    res.trXMinF            = -0.5
    res.vpXF               = 0.20
    res.vpYF               = 0.90

    res.tmYUseLeft         = False

    res.xyMarkerSizeF      = 0.02
    res.xyMarkerThicknessF = [4]
    res.xyMarkers          = [6]

    res.xyLineThicknessF   = [4]

    #-----------------------------------------------------------------------------------------

    res.tiXAxisString      = x_name
    res.tiYAxisFontColor   = 'red'
    res.tiYAxisString      = potential_name

    res.trYMinF            = min(tuple(potential_y_l))
    res.trYMaxF            = max(tuple(potential_y_l))

    res.tmYROn             = False
    res.tmYRMinorOn        = False
    res.tmYRLabelsOn       = False

    res.xyMarkLineMode     = 'Markers'
    res.xyMarkerColor      = 'red'

    plot = Ngl.xy(wks, potential_x_l, potential_y_l, res)
    Ngl.draw(plot)

    res.tiXAxisString      = ''
    res.tiYAxisString      = ''

    res.tmXBOn             = False
    res.tmXBMinorOn        = False
    res.tmXBLabelsOn       = False
    res.tmXTOn             = False
    res.tmXTMinorOn        = False
    res.tmXTLabelsOn       = False
    res.tmYLOn             = False
    res.tmYLMinorOn        = False
    res.tmYLLabelsOn       = False
    res.tmYROn             = False
    res.tmYRMinorOn        = False
    res.tmYRLabelsOn       = False

    res.xyMarkLineMode     = 'Lines'
    res.xyLineColor        = 'red'

    new_x_l = []
    new_y_l = []

    x = potential_x_l[0]
    while x <= potential_x_l[-1]:
          new_x_l.append(x)
          y = polyval(potential_coefficients, [x])
          new_y_l.append(y[0])
          x = x + 0.001

    plot = Ngl.xy(wks, new_x_l, new_y_l, res)
    Ngl.draw(plot)

    #-----------------------------------------------------------------------------------------

    res.tiYAxisSide        = 'Right'

    res.tiXAxisString      = ''
    res.tiYAxisFontColor   = 'blue'
    res.tiYAxisString      = property_name

    res.trYMinF            = min(tuple(property_y_l))
    res.trYMaxF            = max(tuple(property_y_l))

    res.tmXBOn             = False
    res.tmXBMinorOn        = False
    res.tmXBLabelsOn       = False
    res.tmXTOn             = False
    res.tmXTMinorOn        = False
    res.tmXTLabelsOn       = False
    res.tmYLOn             = False
    res.tmYLMinorOn        = False
    res.tmYLLabelsOn       = False
    res.tmYROn             = True
    res.tmYRMinorOn        = True
    res.tmYRLabelsOn       = True

    res.xyMarkLineMode     = 'Markers'
    res.xyMarkerColor      = 'blue'

    plot = Ngl.xy(wks, property_x_l, property_y_l, res)
    Ngl.draw(plot)

    res.tiXAxisString      = ''
    res.tiYAxisString      = ''

    res.tmYROn             = False
    res.tmYRMinorOn        = False
    res.tmYRLabelsOn       = False

    res.xyMarkLineMode     = 'Lines'
    res.xyLineColor        = 'blue'

    new_x_l = []
    new_y_l = []

    x = property_x_l[0]
    while x <= property_x_l[-1]:
          new_x_l.append(x)
          y = polyval(property_coefficients, [x])
          new_y_l.append(y[0])
          x = x + 0.001

    plot = Ngl.xy(wks, new_x_l, new_y_l, res)
    Ngl.draw(plot)

    #-----------------------------------------------------------------------------------------

    res.trYMinF            = min(tuple(potential_y_l))
    res.trYMaxF            = max(tuple(potential_y_l))

    res.xyLineColor        = 'black'
    res.xyLineThicknessF   = [2]

    for i in range(len(psi_squared_save)):
        for j in range(len(psi_squared_save[i])):
            psi_squared_save[i][j] = psi_squared_save[i][j] + energy_save[i]

    res.xyDashPatterns     = [2]
    plot = Ngl.xy(wks, q, psi_squared_save[0], res)
    Ngl.draw(plot)

    res.xyDashPatterns     = [16]
    plot = Ngl.xy(wks, q, psi_squared_save[1], res)
    Ngl.draw(plot)

    #-----------------------------------------------------------------------------------------

    Ngl.frame(wks)
