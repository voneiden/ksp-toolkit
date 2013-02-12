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
import time
import scipy.optimize as so
from pylab import figure, fignum_exists, Circle, linspace, pi, mgrid, sin, cos, zeros, cross, norm, ma, isnan
from celestial import Kerbin
import toolkit
from mpl_toolkits.mplot3d import Axes3D

class Orbit2D:
    ''' Rendering window that displays 2D orbit on reference body '''
    def __init__(self,ref=Kerbin):
        self.ref = ref
        self.figure = figure(figsize=(5,5))
        self.axis = self.figure.gca(projection="rectilinear")#,aspect='equal')
        self.axis.set_aspect('equal')
        #self.axis.set_autoscale_on(True) 
        self.objects = []
        self.figure.canvas.manager.window.after(1000, self.update)

        self.linspace = linspace(0,2*pi,360)
    def update(self):
        #print "UPDATE"
        if not fignum_exists(self.figure.number):
            return
        self.figure.canvas.manager.window.after(1000, self.update)
        if len(self.objects) == 0:
            self.axis.clear()
            self.text = self.axis.text(0.5,0.5,"No objects tracked")
            self.text.set_horizontalalignment("center")
            self.text.set_verticalalignment("center")
            self.figure.canvas.draw()
            return
        
        if toolkit.db.UT:
            t = time.time() - toolkit.db.UTt + toolkit.db.UT
        else:
            t= time.time()
        #figure(self.figure.number)
        self.axis.clear()
        planet = Circle((0,0),600000,color="blue")
        self.axis.add_artist(planet)
        for celestial in self.objects:
            if celestial.ref != self.ref or celestial.orbit == None:
                continue
            xpos = []
            ypos = []
            
            for ta in self.linspace:
                rv = celestial.orbit.plot2D(ta) 
                xpos.append(rv[0])
                ypos.append(rv[1])
                
            
            plotOrbit = self.axis.plot(xpos,ypos,color="red")
            rv,vv = celestial.orbit.eph2D(t)
            #print rv,vv
            self.scat = (rv,vv)
            plotPos = self.axis.scatter(rv[0],rv[1],color="green")
            #print "Pos",rv
        self.axis.relim()
        self.axis.autoscale_view(True,True,True)
        #xlim = self.axis.get_xlim()
        #ylim = self.axis.get_ylim()
        #minlim = min((xlim[0],ylim[0]))
        #maxlim = max((xlim[1],ylim[1]))
        
        #self.axis.set_ylim((minlim,maxlim))
        #self.axis.set_xlim((minlim,maxlim))
        #if xlim[0] > ylim[0]:
        #    self.axis.set_ylim(xlim)
        #else:
        #    self.axis.set_xlim(ylim)
            
        self.figure.canvas.draw()
            #self.axis.plot([random.randint(1,6),random.randint(1,6),random.randint(1,6),random.randint(1,6)])
            
        #draw()
        
            
            
    def track(self,target):
        self.objects.append(target)
        
    def untrack(self,target):
        if target in self.objects:
            self.objects.remove(target)
        
class Orbit3D:
    ''' Displays orbits plot in 3d '''
    def __init__(self,ref=Kerbin):
        self.ref = ref
        self.figure = figure(figsize=(5,5))
        self.axis = self.figure.gca(projection='3d')
        self.axis.set_aspect('equal')
        self.axis.set_axisbelow(True)
        self.figure.canvas.manager.window.after(1000, self.update)
        self.objects = {}
        
        self.linspace = linspace(0,2*pi,360)
        self.sphere = None
        self.notrack = False
    
    def update(self):
        #print "UPDATE"
        if not fignum_exists(self.figure.number):
            return
        self.figure.canvas.manager.window.after(1000, self.update)
        if len(self.objects) == 0:
            self.axis.clear()
            self.text = self.axis.text(0.5,0.5,0.5,"No objects tracked")
            self.text.set_horizontalalignment("center")
            self.text.set_verticalalignment("center")
            self.figure.canvas.draw()
            self.notrack = True
            return
           
        if toolkit.db.UT:
            t = time.time() - toolkit.db.UTt + toolkit.db.UT
        else:
            t= time.time()
        #figure(self.figure.number)
        #self.axis.clear()
        if self.notrack:
            self.notrack = False
            self.axis.clear()
            scale=600000
            u, v = mgrid[0:2*pi:20j, 0:pi:10j]
            x=cos(u)*sin(v)*scale
            y=sin(u)*sin(v)*scale
            z=cos(v)*scale
            self.sphere = self.axis.plot_wireframe(x, y, z, color="r")

        for celestial,plots in self.objects.items():
            if celestial.ref != self.ref or celestial.orbit == None:
                continue
            xpos = []
            ypos = []
            zpos = []
            
            for ta in self.linspace:
                rv = celestial.orbit.plot(ta) 
                
                xpos.append(rv[0])
                ypos.append(rv[1])
                zpos.append(rv[2])
            
            if not plots[0]:
                plots[0] = self.axis.plot(xpos,ypos,zpos,color="red")[0]
            else:
                plots[0].set_data(xpos,ypos)
                plots[0].set_3d_properties(zpos)
                
            
            
            rv,vv = celestial.eph(t)
            #print rv,vv
            #self.scat = (rv,vv)
            if not plots[1]:
                plots[1] = self.axis.scatter(rv[0],rv[1],rv[2],color="green")
            else:
                #plots[1].set_data(rv[0],rv[1],rv[2])
                plots[1].set_offsets([[rv[0],rv[1]]])
                plots[1].set_3d_properties([rv[2]],"z")
                
            #print "Pos",rv
        self.axis.relim()
        
        xlim = self.axis.get_xlim()
        ylim = self.axis.get_ylim()
        
        m = max(abs(xlim[0]),abs(xlim[1]),abs(ylim[0]),abs(ylim[1]))
        lim = (-m,m)
        
        self.axis.set_xlim(lim)
        self.axis.set_ylim(lim)
        self.axis.set_zlim(lim)
        
        #self.axis.autoscale_view(True,True,True,True)
        #xlim = self.axis.get_xlim()
        #ylim = self.axis.get_ylim()
        #minlim = min((xlim[0],ylim[0]))
        #maxlim = max((xlim[1],ylim[1]))
        
        #self.axis.set_ylim((minlim,maxlim))
        #self.axis.set_xlim((minlim,maxlim))
        #if xlim[0] > ylim[0]:
        #    self.axis.set_ylim(xlim)
        #else:
        #    self.axis.set_xlim(ylim)
            
        self.figure.canvas.draw()
            #self.axis.plot([random.randint(
            
    def track(self,target):
        #self.objects.append(target)
        self.objects[target] = [None,None]
        
    def untrack(self,target):
        if target in self.objects:
            del self.objects[target]
        #if target in self.objects:
        #    self.objects.remove(target)
class DistancePlot:
    def __init__(self,a,b):
        if a.ref != b.ref:
            raise AttributeError("Coordinate conversion not implemented yet")
        if not a.isFlying or not b.isFlying:
            raise AttributeError("Both vessels must be orbiting.")
        
        
        self.figure = figure(figsize=(8,8))
        self.axis = self.figure.gca(projection="rectilinear")#,aspect='equal')
        
        
        t = time.time() # Temp
        ts = a.orbit.synodicPeriod(b.orbit)
        
        boundary = [t,t+ts]
        
        X = linspace(boundary[0],boundary[1],1000) #arange(boundary[0],boundary[1],100)
        Y = []
        for i in X:
            Y.append(a.distanceTo(i,b))
        
        self.axis.plot(X,Y)
        d=lambda x: a.distanceTo(x,b)
        closest = so.minimize_scalar(d,method="bounded",bounds=boundary)
        
        self.axis.scatter(closest.x,closest.fun)
        if closest.fun < b.SoI:
            self.axis.plot(boundary,[b.SoI,b.SoI],label="%s SoI"%b.name)
            
            d=lambda x: a.distanceTo(x,b) - b.SoI
            encounter = so.brentq(d,t,closest.x)
            escape = so.brentq(d,closest.x,t+ts)
            
            self.axis.scatter(encounter,b.SoI)
            self.axis.scatter(escape,b.SoI)




class PorkchopPlot:
    def __init__(self, departurePlanet, arrivalPlanet, departureTimes, arrivalTimes, multiplier=86400):
        ''' 
            Generates a porkchop plot
            
            departurePlanet - reference to departure planet
            arrivalPlanet - reference to arrival planet
            departureTimes -  list of departure times
            arrivalTimes - list of arrival times
            multiplier - default 86400 = times are in earth days, set 1 for seconds
        '''
        if departurePlanet.ref != arrivalPlanet.ref:
            print "Lambert solver can only be used on objects with same reference"
            return
        
        mu = departurePlanet.ref.mu
        
        self.figure = figure()
        self.axis = self.figure.gca()
        
        Z = zeros((len(arrivalTimes),len(departureTimes)))
        for depi,deptime in enumerate(departureTimes):
            depsecs = deptime * multiplier
            dp,dv = departurePlanet.eph(depsecs)

            for arri, arrtime in enumerate(arrivalTimes):
                arrsecs = arrtime * multiplier
                ap,av = arrivalPlanet.eph(arrsecs)
                
                dt = arrsecs-depsecs
                # We always want prograde trajectories so..
                if cross(dp,ap)[2] < 0:
                    dt = -dt
                # Solve lambert
                lv1,lv2 = toolkit.lambert(dp,ap,dt,0,mu)
                
                c3 =  norm(lv1-dv) + norm(lv2-av)
                if isnan(c3) or c3 > 2000:
                    c3 = 0
                Z[arri][depi] = c3
                print "v",Z[arri][depi] 
                
        Z = ma.masked_equal(Z,0) 
        self.contour = self.axis.contour(departureTimes,arrivalTimes,Z,colors="k")
        self.contourf = self.axis.contourf(departureTimes,arrivalTimes,Z)
        self.colorbar = self.figure.colorbar(self.contourf)
        
        self.axis.set_ylabel("Arrival day")
        self.axis.set_xlabel("Departure day")
        
