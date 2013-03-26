# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 13:06:09 2013

@author: snaipperi
"""

import  toolkit as tk
from pylab import degrees,arccos, norm
import scipy.optimize as so

def getMohoEvePA(t):
    moho = tk.Moho.eph(t)[0][:2]
    eve = tk.Eve.eph(t)[0][:2]
    
    moho = moho / norm(moho)
    eve = eve / norm(eve)
    
    return degrees(arccos(moho.dot(eve)))
    
    

def getMohoKerbinPA(t):
    moho = tk.Moho.eph(t)[0][:2]
    kerbin = tk.Kerbin.eph(t)[0][:2]
    
    moho = moho / norm(moho)
    kerbin = kerbin / norm(kerbin)
    
    return degrees(arccos(moho.dot(kerbin)))
    
    
    
def getPAdiff(t):
    return abs(getMohoEvePA(t)) + abs(getMohoKerbinPA(t))
    
    
def optimize(sy,ey):
    boundary = [sy*365*24*60*60,ey*365*24*60*60]
    
    match = so.minimize_scalar(getPAdiff,method="bounded",bounds=boundary)
    
    print "Match found at",match.x
    print "PAdiff is",getPAdiff(match.x),"degrees"

def optimizei2(sy,ey):
    boundary = [sy*365*24,ey*365*24]
    
    smallest = 99999999999
    at = 0
    
    for i in xrange(boundary[0],boundary[1]):
        diff = getPAdiff(i*60*60)
        if diff < smallest:
            smallest = diff
            at = i*60*60
            
            
    print "Smallest",smallest
    print "At",at