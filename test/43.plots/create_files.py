#!/usr/bin/env python

import sys, string, os, re

#-------------------------------------------------------------------------------

def to_file(s, file):

    f = open(file.lower(), 'w')
    f.write(s)
    f.close()

#-------------------------------------------------------------------------------

def from_file(file):

    f = open(file, 'r')
    s = f.read()
    f.close()

    return s

#-------------------------------------------------------------------------------

operator_d = {}

operator_d['1x'] = """ 'operator'
 XAVECTOR
 'NEF 03'
 'NEF 02'
 COMFACTOR
 -137.03599980944"""
operator_d['1y'] = """ 'operator'
 YAVECTOR
 'NEF 01'
 'NEF 03'
 COMFACTOR
 -137.03599980944"""
operator_d['1z'] = """ 'operator'
 ZAVECTOR
 'NEF 02'
 'NEF 01'
 COMFACTOR
 -137.03599980944"""

operator_d['2x'] = """ 'operator'
 XAVECTOR
 'NEF 06'
 'NEF 05'
 COMFACTOR
 -137.03599980944"""
operator_d['2y'] = """ 'operator'
 YAVECTOR
 'NEF 04'
 'NEF 06'
 COMFACTOR
 -137.03599980944"""
operator_d['2z'] = """ 'operator'
 ZAVECTOR
 'NEF 05'
 'NEF 04'
 COMFACTOR
 -137.03599980944"""

operator_d['ex'] = " 'XDIPLEN'"
operator_d['ey'] = " 'YDIPLEN'"
operator_d['ez'] = " 'ZDIPLEN'"

operator_d['bx'] = """ 'operator'
 XAVECTOR
 'ZDIPLEN'
 'YDIPLEN'
 COMFACTOR
 -68.517999904721"""
operator_d['by'] = """ 'operator'
 YAVECTOR
 'XDIPLEN'
 'ZDIPLEN'
 COMFACTOR
 -68.517999904721"""
operator_d['bz'] = """ 'operator'
 ZAVECTOR
 'YDIPLEN'
 'XDIPLEN'
 COMFACTOR
 -68.517999904721"""

#-------------------------------------------------------------------------------

template_response = from_file('template_response')

for operator in operator_d:

    s = string.replace(template_response, 'operator', operator_d[operator])

    file = 'response_' + operator + '.inp'
    to_file(s, file)
