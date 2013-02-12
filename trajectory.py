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
import celestial
from pylab import figure, linspace
import scipy.optimize as so

class Trajectory:
    ''' This class calculates trajectory through patched conics.. hopefully!'''
    def __init__(self,epoch,ref,pos,vel,maxdepth=5,steps=200):
        self.epoch = epoch
        self.ref = ref
        self.position = pos
        self.velocity = vel
        self.maxdepth = maxdepth
        self.depth = 0
        self.steps = steps
        
        self.trajectory = []
        
        self.enter(ref,epoch,pos,vel)
        print "Trajectory calculation completed"
        
    def enter(self,ref,epoch,pos,vel):
        print
        print "#"*30
        print "ENTRY POINT (T=%i)"%epoch
        print "#"*30
        
        ''' Begin keplerian simulation at reference body at defined velocities '''
        self.depth += 1
        if self.depth == self.maxdepth:
            print "Max patched conics iterations"
            return
        print "Epoch pos,vel",pos,vel
        orbit = celestial.Orbit(None,ref,epoch,pos,vel)
        print "Entry-point from orbit",orbit.eph3D(orbit.eph2D(epoch))
        period = orbit.period()
        if orbit.e > 1:
            # What shall we use as a reference time when shit is hyperbolic
            period = abs(orbit.period() * 2)
            print "Hyperbolic period (x2):",period
            
        if self.depth > 1:
            boundary = [epoch+5,epoch+period]
        else:
            boundary = [epoch,epoch+period]
            
        
        
        # Things we want to look for
        # 2) Escaping reference body SoI
        # 1) Entering child body SoI
        
        for child in ref.childs:
            d2child = lambda x: orbit.distanceTo(x,child)
            fig = figure()
            fig.canvas.manager.set_window_title('(Depth: %i | Encounter) Distance to %s'%(self.depth,child.name) )
            ax = fig.gca(projection="rectilinear")
            X = linspace(boundary[0],boundary[1],1000)
            Y = []
            for i in X:
                Y.append(d2child(i))
            ax.plot(X,Y)
            ax.plot(boundary,[child.SoI,child.SoI],label="%s SoI"%child.name)
            
            periapsis = so.minimize_scalar(d2child,method="bounded",bounds=boundary)
            print "Boundary",boundary
            print "Min found at",periapsis.x
            if periapsis.fun < child.SoI:
                
                #ax.scatter(periapsis.x,d2child(periapsis.x))
                    
                try: 
                    d2child = lambda x: orbit.distanceTo(x,child) - child.SoI
                    encounter = so.brentq(d2child,boundary[0],periapsis.x)
                   
                except ValueError:
                    print "Did not encounter",child.name
                    pass # Not encountering child SoI
                else:
                    print "Encounter!",encounter
                    ax.scatter(encounter,orbit.distanceTo(encounter,child))
                    self.trajectory.append([orbit,epoch,encounter])
                    pos,vel = orbit.eph3D(orbit.eph2D(encounter))
                    rpos,rvel = celestial.Convert(encounter,ref,child,(pos,vel))
                    print "Converting from",pos,vel
                    print "To (entry-point)",rpos,rvel
                    print "#"*30
                    print "ENCOUNTER POINT (T=%i)"%encounter
                    print "#"*30
                    self.enter(child,encounter,rpos,rvel)
                    return
        
        d2ref = lambda x: orbit.distanceTo(x,ref) - ref.SoI
        if ref.ref != None:
            fig = figure()
            fig.canvas.manager.set_window_title('(Depth: %i | Escape) Distance to %s'%(self.depth,ref.name) )
            ax = fig.gca(projection="rectilinear")
            X = linspace(boundary[0],boundary[1],1000)
            Y = []
            for i in X:
                Y.append(orbit.distanceTo(i,ref))
            ax.plot(X,Y)
            ax.plot(boundary,[ref.SoI,ref.SoI],label="%s SoI"%ref.name)
            print "Boundary",boundary
            try: 
                escape = so.brentq(d2ref,boundary[0],boundary[1])
            except ValueError:
                pass # Not escaping reference SoI
                print "No escape from",orbit.ref.name
            else:
                print "Escape!"
                ax.scatter(escape,orbit.distanceTo(escape,ref))
                self.trajectory.append([orbit,epoch,escape])
                pos,vel = orbit.eph3D(orbit.eph2D(escape))
                rpos,rvel = celestial.Convert(escape,ref,ref.ref,(pos,vel))
                print "#"*30
                print "ESCAPE POINT (T=%i)"%escape
                print "#"*30
                self.enter(ref.ref,escape,rpos,rvel)
                return
                
        self.trajectory.append([orbit,epoch,epoch+period])
        
        
