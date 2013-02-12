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
from math import sqrt
from pylab import array, radians, degrees, cross, linalg, sin, cos, arccos, pi, arctan2
from pylab import arccosh, sinh, cosh, arctan, tanh, arcsin, arcsinh
import scipy.optimize as so

PI2 = pi*2

class Celestial(object):
    def __init__(self,name,ref=None,**kwargs):
        print name,ref,kwargs
        '''
        name - String - Name of the object
        ref - KeplerObject - reference body (None for Kerbol)
        pid - String - Unique ID
        coords - List - [Longitude,Latitude], sets isFlying to False
        elements - List - (ref,epoch,a,e,i,lan,aop,M0), see Orbit class, isFlying True
        isFlying - Bool - override
        mu - Float - standard gravitational parameter of object
        '''
        self.name = name
        self.ref  = ref
        self.childs = []
        self.paths = {}
        self.depth = 0 # For coordinate system conversion
        self.pid  = None
        self.coords = None
        self.orbit = None
        self.isFlying = None
        self.mu = None
        self.radius = None
        self.SoI = None
        keys = kwargs.keys()
        if self.ref != None:
            self.depth = self.ref.depth + 1
            if "coords" not in keys and "elements" not in keys and "state" not in keys:
                raise ValueError("Either coords or elements argument must be specified")
            else:
                
                if "coords" in keys:
                    self.isFlying = False
                    self.coords = kwargs["coords"]
                
                if "elements" in keys:
                    self.isFlying = True
                    self.orbit = Orbit(self,self.ref,*kwargs["elements"])
                if "state" in keys:
                    self.isFlying = True
                    self.orbit = Orbit(self,self.ref,*kwargs["state"])
                    
                if "isFlying" in keys:
                    self.isFlying = kwargs["isFlying"]
            
        if "pid" in keys:
            self.pid = kwargs["pid"]
        if "mu" in keys:
            self.mu = float(kwargs["mu"])
        if "depth" in keys:
            self.depth = kwargs["depth"] 
        if "radius" in keys:
            self.radius = kwargs["radius"] 
        if "SoI" in keys:
            self.SoI = kwargs["SoI"]

        
        print("New Keplerian object created")
        print("- "+self.name)
        if self.ref:
            print("- Orbital period (min): %i"%(self.orbit.period()/60))
    
        
    def eph(self,epoch):
        if self.ref == None:
            return ([0,0,0],[0,0,0])
            
        rv,vv = self.orbit.eph3D(self.orbit.eph2D(epoch))
        
        return (rv,vv)

    def distanceTo(self,epoch,obj):
        return self.orbit.distanceTo(epoch,obj)

    def generatePaths(self):
        ''' This function generates path to every object in the parent-child tree
            the information is used in coordinate conversion. This function should
            be called once all major celestial bodies have been generated/updated.'''
        
        self.paths = {}
        closedNodes = []
        openNodes = self.childs[:] + [self.ref]
        
        while len(openNodes) > 0:
            node = openNodes.pop()
            if node == None:
                continue
            
            if node.ref == self: # Case 1: Node is self child
                self.paths[node] = [node]
            elif self in node.childs: # Case 2: Node is self parent
                self.paths[node] = [node]
            elif node.ref in self.paths.keys():# Case 3: Node is child of another known node
                self.paths[node] =  self.paths[node.ref] + [node]
            else:
                for childnode in node.childs: # Case 4: Node is parent of
                    if childnode in self.paths.keys():
                        self.paths[node] = self.paths[childnode] + [node]
                        break
            if node not in self.paths[node]:
                print "Error"
                print node,node.name
                raise RuntimeError("Unable to determine node")
            
            closedNodes.append(node)
            for newnode in node.childs + [node.ref]:
                if newnode != None and newnode not in closedNodes and newnode not in openNodes and newnode != self:
                    openNodes.append(newnode)        
            
    
class Coordinate:
    def __init__(self,ref,position,velocity):
        self.position = position
        self.velocity = velocity
        self.ref = ref
        
    def __str__(self):
        print "<Reference: %s - Position:%s>"%(self.ref.name,str(self.coordinate))
        
        
def Convert(epoch,ref,toref,components):
    ''' Converts components from being relative to fref
        to being relative to tref at given time t'''       
    
    
    if toref == ref:
        return components
        
    lastdepth = ref.depth
    lastref   = ref
    position = array(components[0])
    velocity = array(components[1])
    print "Convert()ing from reference",ref.name
    print "To referece",toref.name
    print ref.paths.keys()
    for celestial in ref.paths[toref]:
        if celestial.depth < lastdepth:
            lastref_r,lastref_v = lastref.eph(epoch)
            print "Moving towards parent",lastref_r,lastref_v
            #position = [position[0] + lastref_r[0], position[1] + lastref_r[1], position[2] + lastref_r[2]]
            #velocity = [velocity[0] + lastref_v[0], velocity[1] + lastref_v[1], velocity[2] + lastref_v[2]]
            position[0] = position[0] + lastref_r[0]
            position[1] = position[1] + lastref_r[1]
            position[2] = position[2] + lastref_r[2]
            velocity[0] = velocity[0] + lastref_v[0]
            velocity[1] = velocity[1] + lastref_v[1]
            velocity[2] = velocity[2] + lastref_v[2]
        
        elif celestial.depth > lastdepth:
            celestial_r,celestial_v = celestial.eph(epoch)
            print "Moving towards child",celestial_r,celestial_v
            #position = [position[0] - celestial_r[0], position[1] - celestial_r[1], position[2] - celestial_r[2]]
            #velocity = [velocity[0] - celestial_v[0], velocity[1] - celestial_v[1], velocity[2] - celestial_v[2]]
            position[0] = position[0] - celestial_r[0]
            position[1] = position[1] - celestial_r[1]
            position[2] = position[2] - celestial_r[2]
            velocity[0] = velocity[0] - celestial_v[0]
            velocity[1] = velocity[1] - celestial_v[1]
            velocity[2] = velocity[2] - celestial_v[2]
            
        lastdepth = celestial.depth
        lastref = celestial
        
    return (position,velocity)
    
            
            
class Orbit(object):
    def __init__(self,obj,ref,epoch,*args):
        '''
        *args are either a,e,i,lan,aop,M0
        or r,v
        '''
        self.obj = obj
        self.ref = ref
        self.shape2D = None
        self.shape3D = None
        self._a = None
        self._e = None
        
        if len(args) == 2:
            self.initFromStateVectors(epoch,*args)
        elif len(args) == 6:
            self.epoch = float(epoch)
            self._a = float(args[0])
            self._e = float(args[1])
            self.i = radians(float(args[2]))
            self.lan = radians(float(args[3]))
            self.aop = radians(float(args[4]))
            self.M0 = float(args[5])
            self.calculateH()
            #if self.e < 1:
            #    self.h = sqrt(self.a*self.ref.mu*(1-self.e**2))
            #elif self.e > 1:
            #    self.h = sqrt(-self.a*self.ref.mu*(self.e**2-1))
        else:
            raise AttributeError("Unable to init"+str(args))
        if len(args) == 0:
            pass
        else:
            self.updateShape()
    
    @property
    def a(self): return self._a
    
    @a.setter
    def a(self,value):
        self._a = value
        if self.e != None and self.a != None:
            self.calculateH()
    
    @property
    def e(self): return self._e
    
    @e.setter
    def e(self,value):
        self._e = value
        if self.e != None and self.a != None:
            self.calculateH()
        
    def calculateH(self):
        print "RECALCULATING H"
        if self.e < 1 and self.a > 0:
            self.h = sqrt(self.a*self.ref.mu*(1-self.e**2))
        elif self.e > 1 and self.a < 0:
            self.h = sqrt(-self.a*self.ref.mu*(self.e**2-1))
            
        elif self.e == 0:
            raise AttributeError("Parabolic trajectory, unable to proceed")

        else:
            self.h = None
    def initFromStateVectors(self,epoch,pV,vV):
        self.epoch = epoch
        
        # 1) Calculate auxilary vector h
        hV = cross(pV,vV)
        
        
        # 2) Normalize position,velocity, specific angular momentum, calculate radial velocity 
        
        p = linalg.norm(pV)
        v = linalg.norm(vV)
        h = linalg.norm(hV)
        print "H:",h
        radv = pV.dot(vV) / p
        hVu = hV / h
        pVu = pV / p
        nV = cross(array([0,0,1]),hV)
        n = linalg.norm(nV)
        if n == 0:
            nVu = array([0,0,0])
        else:
            nVu = nV/n
        # 3) Calculate inclination
        #self.i = arccos(hV[2]/h)
        self.i = arcsin(linalg.norm(cross(array([0,0,1]),hVu)))
        print "i1",self.i
        print "RADVEL",radv
        self.i = arccos(array([0,0,1]).dot(hV)/h)
        #if radv < 0:
        #    self.i = PI2 - self.i 
        print "i2",self.i
        # 4) Calculate node line
        
        
        # 5) Calculate longitude of ascending node = right ascension of ascending node
        '''
        if self.i == 0:
            self.lan=0
        elif nV[1] >= 0:
            self.lan = arccos(nV[0] / n)
        else:
            self.lan = PI2 - arccos(nV[0] / n)
        '''
        
        if self.i == 0:
            self.lan = 0
        else:
            self.lan = arcsin(cross(array([1,0,0]),nVu).dot(array([0,0,1])))
            print "lan1",self.lan
            self.lan = arccos(array([1,0,0]).dot(nV)/n)
            if nV[1] < 0:
                self.lan = PI2-self.lan
            print "lan2",self.lan
        
        # 6) Eccentricity vector
        #eV = (1.0 / self.ref.mu)*((v**2 - (self.ref.mu / p))*pV - radv*vV)
        #eV2 = (1.0 / self.ref.mu) * ( hV - self.ref.mu * (pV/p))
        #eV3 = hV/self.ref.mu - (pV/p)
        
        # Source: cdeagle
        eV = cross(vV,hV)/self.ref.mu - pVu
        #print "eV1:",eV,linalg.norm(eV)
        #print "eV2:",eV2,linalg.norm(eV2)
        #print "eV3:",eV3,linalg.norm(eV3)
        print "eV3:",eV,linalg.norm(eV)
        self._e = linalg.norm(eV)
        #eVu = eV / self.e
        
        print "h",h
        print "u",self.ref.mu
        print "v",v
        print "r",p
        
        print "alte:",sqrt(1+(h**2/self.ref.mu**2)*(v**2-(2*self.ref.mu)/p)**2)
        # 7) Argument of perigree
        '''
        if self.e == 0:
            self.aop = 0
        elif self.i == 0:
          self.aop = arccos(eV[0] / self.e)  
        elif eV[2] >= 0:
            print "AOP AOP AOP"
            #self.aop = arccos(nV.dot(eV) / (n*self.e))
            print cross(nV,eV).dot(hV)
            self.aop = arcsin(cross(nVu,eVu).dot(hVu))
            #self.aop = arccos(n*self.e)
        else:
            self.aop = PI2 - arccos(nV.dot(eV) / (n*self.e))
        '''
        #CDEagle method
        # TODO CHECK how KSP handles this. 
        if self.e == 0:
            self.aop = 0
        elif self.i == 0 and self.e != 0:            
            #self.aop = arccos(eV[0] / self.e)
            #self.aop = arctan2(eV[1],eV[0])
            self.aop = arccos(array([1,0,0]).dot(eV) / self.e)
            print eV
            if eV[2] < 0:
                #self.aop = -self.aop
                self.aop = PI2-self.aop
            
            #print "BOOM",eV
            #if eV[2] < 0:
            #    print "BAM NIGGA"
            #    self.aop = PI2 - self.aop
        elif self.i == 0 and self.e == 0:
            #raise AttributeError("Perfectly circular orbits are not supported atm")
            self.aop = 0
        else:
            #self.aop = arcsin(cross(nVu,eVu).dot(hVu))
            self.aop = arccos(nV.dot(eV)/(n*self.e))
            if eV[2] < 0:
                self.aop = PI2-self.aop
        
        # 8) Semi major axis
        aE = v**2/2.0 - self.ref.mu / p
        self._a = -self.ref.mu / (2*aE)
        print "Old method for semi-major",self.a
        self._a = h**2 / (self.ref.mu * (1-self.e**2))
        print "New method for semi-major",self.a  
        #if self.e > 1:
        #    self._a = h**2 / (self.ref.mu * (self.e**2 - 1))
        
        if self.e == 0:
            if self.i == 0: #TODO update document to this
                print "JEA JEA JEA JEA"*10
                
                ta = arccos(array([1,0,0]).dot(pV) / p)
                if pV[1] < 0: # Vallado pg. 111
                    ta = PI2 - ta
            else: #TODO VERIFY THIS CASE
                ta = arccos((nV.dot(pV))/(n*p))
                if pV[2] < 0: # Vallado pg. 110
                    ta = PI2 - ta
            E = ta
            self.M0 = E
              
        elif self.e < 1:
            # 9) True anomaly, eccentric anomaly and mean anomaly
            if radv >= 0:
                ta = arccos((eV / self.e).dot(pV/p))
            else:
                ta = PI2 - arccos((eV / self.e).dot(pV/p))
            
            
            E = arccos((self.e+cos(ta))/(1+ self.e*cos(ta)))
            if radv < 0:
                E = PI2 - E
        
            self.M0 = E - self.e * sin(E)
            
        elif self.e > 1:
            # 9) Hyperbolic True anomaly, eccentric anomaly and mean anomaly
            # http://scienceworld.wolfram.com/physics/HyperbolicOrbit.html
            V = arccos((abs(self.a)*(self.e**2 - 1)) /(self.e * p) - 1/self.e)
            ta = arccos((eV / self.e).dot(pV/p))
            
            if radv < 0: #TODO: Should affect F too?
                # Negative = heading towards periapsis
                print "PI2"
                V = PI2 - V
                ta = PI2-ta
            print "V",V
            print "TA",ta
            # http://www.bogan.ca/orbits/kepler/orbteqtn.html In you I trust
            # Hyperbolic eccentric anomaly
            cosV = cos(V)
            F = arccosh((self.e+cosV)/(1+self.e*cosV))
            if radv < 0:
                F = -F
            F2 = arcsinh((sqrt(self.e-1)*sin(V))/(1+self.e*cos(V)))
            ##F1 = F2
            print "F1:",F
            print "F2:",F2
            self.M0 = self.e * sinh(F) - F
            
        
     
        
        
        self.h = h
        
        print "Semi-major:",self.a
        print "Eccentricity:",self.e
        print "Inclination:",degrees(self.i),"deg"
        print "LAN:",degrees(self.lan),"deg"
        print "AoP:",degrees(self.aop),"deg"
        print "Mean anomaly:",self.M0
        print "Specific angular momentum:",self.h
        if self.e < 1:
            print "Eccentric anomaly",E
            print "True anomaly",ta
        else:
            print "Hyperbolic eccentric anomaly",F
            print "Hyperbolic true anomaly",degrees(V)
            
        print "Distance from object:",p
        print "Velocity:",v
        
        


        
    def distanceTo(self,epoch,obj):        
        if self.ref == obj:
            op,ov = ([0,0,0],[0,0,0])
        elif self.ref != obj.ref:
            print self.ref.name, "->", obj.ref.name , "(",obj.name,")"
            raise AttributeError("Coordinate conversion not implemented")
        else:
            op,ov = obj.eph(epoch)
        sp,sv = self.eph3D(self.eph2D(epoch))
        
        return linalg.norm(sp-op)    
    

    def updateShape(self):
        ''' This function needs to implement
            patched conics.. '''
        return
        if self.ref == None:
            return
        period = self.period()
        
        steps = 200
        steptime = period / steps
        
        X2D = []
        Y2D = []
        X3D = []
        Y3D = []
        Z3D = []
        #Z = []
        for i in xrange(steps+1):
            
            rv,vv = self.eph2D(steptime*i)
            #print "Got",rv,vv
            X2D.append(rv[0])
            Y2D.append(rv[1])
            
            rv,vv = self.eph3D((rv,vv))
            X3D.append(rv[0])
            Y3D.append(rv[1])
            Z3D.append(rv[2])


        self.shape2D = [X2D,Y2D]
        self.shape3D = [X3D,Y3D,Z3D]
    def period(self):
        ''' Returns the period of current orbit '''
        if self.e < 1:
            return PI2*sqrt(self.a**3/self.ref.mu)
        elif self.e > 1:
            return 2 * self.M0 * sqrt(-self.a**3 / self.ref.mu)
        else:
            raise AttributeError("Eccentricity 0 not defined")
            return None
    
    def synodicPeriod(self,orbit):
        if orbit.ref != self.ref:
            raise AttributeError("Synodic period can only be calculated for orbits with the same reference")
        return 1.0 / abs(1.0/self.period() - 1.0 / orbit.period())
    
    def plot(self,ta):
        return self.eph3D((self.plot2D(ta),None))
        
    def plot2D(self,ta):
        r = (self.h**2/self.ref.mu)*(1.0/(1.0+self.e*cos(ta)))
        rv = r * array([cos(ta),sin(ta),0])
        
        #v = self.ref.mu / self.h
        #vv = v * array([-sin(ta),self.e+cos(ta),0])
        
        return rv#,vv)
    

        
        
    def eph2D(self,epoch):
        
        if self.ref == None:
            return ([0,0,0],[0,0,0])
                
        dt = epoch-self.epoch  
        #print "dT",dt
        
        # Step 1 - Determine mean anomaly at epoch
        if self.e == 0:
            M = self.M0 + dt * sqrt(self.ref.mu / self.a**3)
            M %= PI2
            E = M
            ta = E
            r3 = (self.h**2/self.ref.mu)
            rv = r3 * array([cos(ta),sin(ta),0])
            v3 = self.ref.mu / self.h
            vv = v3 * array([-sin(ta),self.e+cos(ta),0])
            
            return (rv,vv)
            
        if self.e < 1:
            if epoch == self.epoch:
                M = self.M0
            else:
                M = self.M0 + dt * sqrt(self.ref.mu / self.a**3)
                M %= PI2
             
             # Step 2a - Eccentric anomaly
            E = so.newton(lambda x: x-self.e * sin(x) - M,M)
            
            # Step 3a - True anomaly, method 1
            ta = 2 * arctan2(sqrt(1+self.e)*sin(E/2.0), sqrt(1-self.e)*cos(E/2.0))
            #print "Ta:",ta
            # Method b is faster
            cosE = cos(E)
            
            ta2 = arccos((cosE - self.e) / (1-self.e*cosE))
            #print "e",self.e
            #print "M",M
            #print "E",E
            #print "TA:",ta
            #print "T2:",ta2
            
            # Step 4a - distance to central body (eccentric anomaly).
            r = self.a*(1-self.e * cos(E))
            
            # Alternative (true anomaly)
            r2 = (self.a*(1-self.e**2) / (1.0 + self.e * cos(ta)))
            
            # Things get crazy
            #h = sqrt(self.a*self.ref.mu*(1-self.e**2))
            r3 = (self.h**2/self.ref.mu)*(1.0/(1.0+self.e*cos(ta)))
            
        
            #print "R1:",r
            #print "R2:",r2
            #print "R3:",r3
            rv = r3 * array([cos(ta),sin(ta),0])
            
            #v1 = sqrt(self.ref.mu * (2.0/r - 1.0/self.a))
            #v2 = sqrt(self.ref.mu * self.a) / r
            v3 = self.ref.mu / self.h
            #meanmotion = sqrt(self.ref.mu / self.a**3)
            #v2 = sqrt(self.ref.mu * self.a)/r
            
            #print "v1",v1
            #print "v2",v2
            #print "v3",v3
            #print "mm",meanmotion
            
            
            # Velocity can be calculated from the eccentric anomaly            
            #vv = v1 * array([-sin(E),sqrt(1-self.e**2) * cos(E),0])
            
            # Or from the true anomaly (this one has an error..)
            #vv = sqrt(self.ref.mu * self.a)/r * array([-sin(ta),self.e+cos(ta),0])
            
            vv = v3 * array([-sin(ta),self.e+cos(ta),0])
            
            #print "rv",rv
            #print "vv",vv
            #print "check",E,-sin(E),v1*-sin(E)
            #print "V1:",vv
            #print "V2:",vv2

            return (rv,vv)
            
        elif self.e > 1:
            # Hyperbolic orbits
            # Reference: Stephen Kemble: Interplanetary Mission Analysis and Design, Pg 220-221
            M = self.M0 + dt * sqrt(-(self.ref.mu / self.a**3))
            #print "M0",self.M0
            #print "M",M
            #print "M",M
            #print "M0",self.M0
            # Step 2b - Hyperbolic eccentric anomaly
            #print "Hyperbolic mean anomaly:",M
            F = so.newton(lambda x: self.e * sinh(x) - x - M,M,maxiter=1000)
            #F = -F
            H = M / (self.e-1)
            #print "AAAA"*10
            #print "F:",F
            #print "H:",H
            #F=H
            #print "Hyperbolic eccentric anomaly:",F
            
            # Step 3b - Hyperbolic true anomaly?
            # This is wrong prooobably            
            #hta = arccos((cosh(F) - self.e) / (1-self.e*cosh(F)))
            #hta = arccos((self.e-cosh(F)) / (self.e*cosh(F) - 1))
            # TÄSSÄ ON BUGI
            hta = arccos((cosh(F) - self.e) / (1 - self.e*cosh(F)))
            hta2 = 2 * arctan2(sqrt(self.e+1)*sinh(F/2.0),sqrt(self.e-1)*cosh(F/2.0))
            hta3 = 2 * arctan2(sqrt(1+self.e)*sinh(F/2.0),sqrt(self.e-1)*cosh(F/2.0))
            hta4 = 2 * arctan2(sqrt(self.e-1)*cosh(F/2.0),sqrt(1+self.e)*sinh(F/2.0))
            #print "Hyperbolic true anomaly:",degrees(hta)
            # This is wrong too            
            #hta2 = 2 * arctan2(sqrt(1+self.e)*sin(F/2.0), sqrt(1-self.e)*cos(F/2.0))
            if M == self.M0:
                print "HTA1:",degrees(hta)
                print "HTA2:",degrees(hta2)
                print "HTA3:",degrees(hta3)
                print "HTA4:",degrees(hta4)
            
            # According to http://mmae.iit.edu/~mpeet/Classes/MMAE441/Spacecraft/441Lecture17.pdf
            # this is right..
            hta = hta2
            #print cos(hta)
            #print cosh(hta)
            
            #####hta = arctan(sqrt((self.e-1.0)/self.e+1.0) * tanh(F/2.0)) / 2.0
            #print "Mean anomaly:",M
            #print "Hyperbolic eccentric anomaly:",F
            #print "HTA:",hta
            #raise
            # Step 4b - Distance from central body?
            # Can calculate it from hyperbolic true anomaly..
            #p = self.a*(1-self.e**2)
            #r = p / (1+self.e * cos(hta))
            r3 = (self.h**2/self.ref.mu)*(1.0/(1.0+self.e*cos(hta)))
            v3 = self.ref.mu / self.h
            # But it's faster to calculate from eccentric anomaly
            #r = self.a*(1-self.e * cos(F))
            
            #print "Hyper r1:",r
            #print "Hyper r2:",r2
            
            rv = r3 * array([cos(hta),sin(hta),0])
            #http://en.wikipedia.org/wiki/Hyperbola
            #rv = array([ r * sin(hta),-self.a*self.e + r * cos(hta), 0])
            #sinhta = sin(hta)
            #coshta = cos(hta)
            #print self.ref.mu,r,self.a,hta
            #vv = sqrt(self.ref.mu *(2.0/r - 1.0/self.a)) * array([-sin(hta),self.e+cos(hta),0])
            vv = v3 * array([-sin(hta),self.e+cos(hta),0])
            
            return (rv,vv)
            
            #raise AttributeError("Oh snap. Hyperbolic orbits..")
        #print "Mean epoch:",M,sqrt(self.ref.mu / self.a**3),self.ref.mu / self.a 
        #print "a**3",self.a**3
        #print "mu",self.ref.mu
        
        
       
        
    def eph3D(self,components):
        lancos = cos(self.lan)
        lansin = sin(self.lan)
        inccos = cos(self.i)
        incsin = sin(self.i)
        argcos = cos(self.aop)
        argsin = sin(self.aop)
        """
        
        
        LAN = array([[ lancos, lansin,0],
                     [-lansin, lancos,0],
                     [0,       0,     1]])
                        
        INC = array([[1,  0,      0],
                     [0,  inccos, incsin],
                     [0, -incsin, inccos]])
                        
        ARG = array([[ argcos, argsin, 0],
                     [-argsin, argcos, 1],
                     [ 0,      0,      1]])
                        
        ROT = LAN.dot(INC).dot(ARG)
        IROT = liqwnalg.inv(ROT)
        """
        # For faster performance construct the inverse
        # transformation matrix straight away
        
        #ROT = array([[-lansin * inccos * argsin + lancos * argcos, lancos * inccos * argsin + lansin * argcos, incsin * argsin],
        #              [-lansin * inccos * argcos - lancos * argsin, lancos * inccos * argcos - lansin * argsin, incsin * argcos],
        #              [ lansin * incsin,                           -lancos*incsin,                              inccos]])
        
        IROT = array([[-lansin * inccos * argsin + lancos * argcos, -lansin * inccos * argcos - lancos * argsin,  lansin * incsin ],
                      [ lancos * inccos * argsin + lansin * argcos,  lancos * inccos * argcos - lansin * argsin, -lancos * incsin],
                      [ incsin * argsin,                             incsin * argcos,                             inccos]])
                     
        #ROTMAT = array([[
        #]])
        #ROT2 = ROT.transpose()
        #TODO the rotation matrix has something silly in it
        #p = array([[components[0][0]],[components[0][1]],[components[0][2]]])
        #
        if components[1] == None:
            return IROT.dot(components[0])
        else:
            return (IROT.dot(components[0]),IROT.dot(components[1]))


class Sun(Celestial):
    def __init__(self,name,**kwargs):
        Celestial.__init__(self,name,None,**kwargs)

class Planet(Celestial):
    def __init__(self,name,ref,**kwargs):
        Celestial.__init__(self,name,ref,depth=ref.depth+1,**kwargs)
        ref.childs.append(self)
        
class Moon(Planet):
    pass
 
class Ship(Celestial):
    def __init__(self,name,ref,**kwargs):
        if ref:
            Celestial.__init__(self,name,ref,depth=ref.depth+1,**kwargs)
        else:
            Celestial.__init__(self,name,ref,**kwargs)

#Define constants
Kerbol = Sun("Kerbol",mu=1.1723328e18,radius=261600000)

Kerbin = Planet("Kerbin", Kerbol,
                elements=[0,13599840256,0,0,0,0,3.14000010490417],
                mu=3531600000000,
                radius=600000,
                SoI=84159286)
                
Mun = Moon("Mun",Kerbin,
           elements=[0,12000000,0,0,0,0,1.70000004768372], 
           mu=65138398000,
           radius=200000,
           SoI=2429559.1)

Minmus = Moon("Minmus",Kerbin,
              elements=[0,47000000,0,6,78,38,0.899999976158142],
              mu=1765800000 ,
              radius=60000,
              SoI=2247428.4 )               
                
            

Duna = Planet("Duna",Kerbol,
              elements=[0,20726155264,
                        0.0509999990463257,
                        0.0599999986588955,
                        135.5,
                        0,
                        3.14000010490417],
               mu=301363210000.0,
               radius=320000.0,
               SoI=47921949.0)
 
 
Celestials = [Kerbol,Kerbin,Mun,Minmus,Duna]
for celestial in Celestials:
    celestial.generatePaths()
    
def testShip():
    return Ship("Testship",Kerbin,elements=[0,8000000,0.5,0,0,0,0])
    
def testElliptic():
    ship = testShip()
    ship.orbit.e = 0.2
    r,v = ship.eph(15)
    print "EPH OUTPUTS",r,v
    print "ORBITEPH OUTPUTS",ship.orbit.eph2D(15)
    print "ORBIT3D OUTPUTS",ship.orbit.eph3D(ship.orbit.eph2D(15))
    print "-"*10
    shipo = Orbit(None,ship.ref,15,r,v)
    print u"µ",ship.orbit.ref.mu, shipo.ref.mu
    print "a",ship.orbit.a, shipo.a
    print "e",ship.orbit.e, shipo.e
    print "i",ship.orbit.i, shipo.i 
    print "l",ship.orbit.lan, shipo.lan
    print "p",ship.orbit.aop, shipo.aop
    print "m",ship.orbit.M0, shipo.M0 
    print "t",ship.orbit.epoch, shipo.epoch 
    sr,sv = shipo.eph3D(shipo.eph2D(15))
    print "r1",r
    print "r2",sr
    print "v1",v
    print "v2",sv

def testElliptic():
    ship = testShip()
    ship.orbit.e = 0.5
    r,v = ship.eph(30000)
    print "EPH OUTPUTS",r,v
    print "ORBITEPH OUTPUTS",ship.orbit.eph2D(15)
    print "ORBIT3D OUTPUTS",ship.orbit.eph3D(ship.orbit.eph2D(15))
    print "-"*10
    shipo = Orbit(None,ship.ref,15,r,v)
    print u"µ",ship.orbit.ref.mu, shipo.ref.mu
    print "a",ship.orbit.a, shipo.a
    print "e",ship.orbit.e, shipo.e
    print "i",ship.orbit.i, shipo.i 
    print "l",ship.orbit.lan, shipo.lan
    print "p",ship.orbit.aop, shipo.aop
    print "m",ship.orbit.M0, shipo.M0 
    print "t",ship.orbit.epoch, shipo.epoch 
    sr,sv = shipo.eph3D(shipo.eph2D(15))
    print "r1",r
    print "r2",sr
    print "v1",v
    print "v2",sv 
    
def testElliptic2(epoch=10000):
    ship = testShip()
    ship.orbit.e = 0.5
    r,v = ship.eph(epoch)
    print "EPH OUTPUTS",r,v
    epoch = epoch
    print "ORBITEPH OUTPUTS",ship.orbit.eph2D(epoch) #NOTE why does this not affect?
    print "ORBIT3D OUTPUTS",ship.orbit.eph3D(ship.orbit.eph2D(epoch))
    print "-"*10
    shipo = Orbit(None,ship.ref,epoch,r,v)
    print u"µ",ship.orbit.ref.mu, shipo.ref.mu
    print "a",ship.orbit.a, shipo.a
    print "e",ship.orbit.e, shipo.e
    print "i",ship.orbit.i, shipo.i 
    print "l",ship.orbit.lan, shipo.lan
    print "p",ship.orbit.aop, shipo.aop
    print "m",ship.orbit.M0, shipo.M0 
    print "t",ship.orbit.epoch, shipo.epoch 
    sr,sv = shipo.eph3D(shipo.eph2D(epoch))
    print "r1",r
    print "r2",sr
    print "v1",v
    print "v2",sv 
    
def testKerbin(epoch=10000):
    ship = Kerbin
    #ship.orbit.e = 0.5
    r,v = ship.eph(epoch)
    print "EPH OUTPUTS",r,v
    epoch = epoch
    print "ORBITEPH OUTPUTS",ship.orbit.eph2D(epoch) #NOTE why does this not affect?
    print "ORBIT3D OUTPUTS",ship.orbit.eph3D(ship.orbit.eph2D(epoch))
    print "-"*10
    shipo = Orbit(None,ship.ref,epoch,r,v)
    print u"µ",ship.orbit.ref.mu, shipo.ref.mu
    print "a",ship.orbit.a, shipo.a
    print "e",ship.orbit.e, shipo.e
    print "i",ship.orbit.i, shipo.i 
    print "l",ship.orbit.lan, shipo.lan
    print "p",ship.orbit.aop, shipo.aop
    print "m",ship.orbit.M0, shipo.M0 
    print "t",ship.orbit.epoch, shipo.epoch 
    sr,sv = shipo.eph3D(shipo.eph2D(epoch))
    print "r1",r
    print "r2",sr
    print "v1",v
    print "v2",sv 
    ship2 = Ship("Foo",None)
    ship2.ref = Kerbol
    ship2.orbit = shipo
    return ship,ship2
    
def testHyperbolic():
    ship = testShip()
    ship.orbit.e = 1.5
    ship.orbit.a = -ship.orbit.a
    r,v = ship.eph(15)
    print "EPH OUTPUTS",r,v
    print "ORBITEPH OUTPUTS",ship.orbit.eph2D(15)
    print "ORBIT3D OUTPUTS",ship.orbit.eph3D(ship.orbit.eph2D(15))
    print "-"*10
    shipo = Orbit(None,ship.ref,15,r,v)
    print u"µ",ship.orbit.ref.mu, shipo.ref.mu
    print "a",ship.orbit.a, shipo.a
    print "e",ship.orbit.e, shipo.e
    print "i",ship.orbit.i, shipo.i 
    print "l",ship.orbit.lan, shipo.lan
    print "p",ship.orbit.aop, shipo.aop
    print "m",ship.orbit.M0, shipo.M0 
    print "t",ship.orbit.epoch, shipo.epoch 
    sr,sv = shipo.eph3D(shipo.eph2D(15))
    print "r1",r
    print "r2",sr
    print "v1",v
    print "v2",sv
    
    
def testEarth():
    earth = Sun("Earth",mu=398600,radius=6000)
    p = array([-4040.0,4815.0,3629.0])
    v = array([-10.39,-4.772,1.744])
    
   #eship1 = Ship("Testship",earth,elements=[0,-16725.20488375983,1.4,0,0,0,0])    
    
    eship = Ship("Testship",earth,state=[0,p,v])
    #eship.orbit.aop = radians(60)
    print eship.eph(0)
    return eship
    
def testMun():
    p = array([ 2413828.03900541,  -276028.29294129,        0.        ])
    v = array([-186.88780283,  241.00746053,    0.00       ])
    
    eship = Ship("Testship",Mun,state=[0,p,v])
    print eship.eph(0)
    print [0,p,v]
    return eship
