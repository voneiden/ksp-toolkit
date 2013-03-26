#!/usr/bin/python2
# -*- coding: utf-8 -*-

"""
   KSP Mission Control
   Copyright (C) 2013 Matti Eiden
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# Yes, this is very evil thing to do in python!
import matplotlib
matplotlib.rcParams['toolbar'] = "None"
import celestial
import trajectory
import plotter
import network
import lambert

# Reload is done to ensure that any changes to the modules are recorded
reload(celestial)
reload(trajectory)
reload(plotter)
reload(network)
reload(lambert)
# And then the namespace is polluted
from celestial import *
from trajectory import *
from plotter import *
from network import *
from lambert import *

class Database(object):
    ''' The database holds information about dynamic content '''
    def __init__(self):
        self.active = None #Active vessel
        self.vessels = {}
        self.refs = {
            "Kerbin":Kerbin,
            "Mun":Mun,
            "Minmus":Minmus,
            "Duna":Duna}
        self.UT = None
        self.UTt = None

db = Database()



print("Toolkit loaded.")
