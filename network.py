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
import socket
import pylab as pl
import toolkit, time
import Tkinter

class Network:
    def __init__(self,ip="192.168.1.13",port=11000):
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        
        #self.figure = pl.figure(figsize=(1,2))
        #self.axis = self.figure.gca(projection="rectilinear")
        self.toplevel = Tkinter.Toplevel()
        self.toplevel.title("Network")
        self.statusmsg = Tkinter.StringVar()
        self.statusmsg.set("Initializing")
        self.label = Tkinter.Label(self.toplevel, textvariable=self.statusmsg)
        self.label.pack()
        
        self.buf = ""
        try:        
            self.socket.connect((ip,port))
            
        except:
            #raise
            #self.text = self.axis.text(0.5,0.5,"Network: Failed")
            #self.text.set_horizontalalignment("center")
            #self.text.set_verticalalignment("center")
            #self.figure.canvas.draw()
            self.statusmsg.set("Connection failed")
        else:
            self.socket.setblocking(0)
            #self.text = self.axis.text(0.5,0.5,"Network: Connected")
            #self.text.set_horizontalalignment("center")
            #self.text.set_verticalalignment("center")
            
            #self.figure.canvas.manager.window.after(1000, self.update)
            #self.figure.canvas.draw()
            self.statusmsg.set("Connected")
            self.after = self.toplevel.after(1000,self.update)
            #print("FUUU")
            
    def update(self):
        #print ("UPDATE NETWORK")
        try:
            data = self.socket.recv(1024)
        except socket.error:
            #print "SocketError"
            # TODO: handle disconnects
            pass
        else:
            if not data:
                self.statusmsg.set("Disconnected")
                #self.text = self.axis.text(0.5,0.5,"Disconnected")
                #self.text.set_horizontalalignment("center")
                #self.text.set_verticalalignment("center")
                #self.figure.canvas.draw()
                #return
            else:
                print "Processing data"
                self.buf += data
                lines = self.buf.split('\n')
                self.buf = lines[-1]
                for line in lines[:-1]:
                    tok = line.split('\t')
                    if tok[0] == "UT":
                        #print "Universal time:",float(tok[1])
                        UT = float(tok[1])
                        toolkit.db.UT = UT
                        toolkit.db.UTt = time.time()
                    
                    elif tok[0] == "V":
                        print "VESSEL", tok[1] 
                        if tok[1] == "D":
                            pid = tok[2]
                            name = tok[3]
                            print "Define vessel (%s) %s"%(tok[2],tok[3])
                            toolkit.db.vessels[pid] = toolkit.Ship(name,None)
                            
                        elif tok[1] == "L":
                            pid = tok[2]
                            ship = toolkit.db.vessels[pid]
                            toolkit.db.UT = float(tok[3])
                            toolkit.db.UTt = time.time() 
                            
                            ship.info["missionTime"]            = float(tok[4])
                            ship.info["acceleration"]           = float(tok[5])
                            ship.info["altitude"]               = float(tok[6])
                            ship.info["angularMomentum"]        = float(tok[7])
                            ship.info["angularVelocity"]        = float(tok[8])
                            ship.info["atmDensity"]             = float(tok[9])
                            ship.info["geeForce"]               = float(tok[10])
                            ship.info["geeForce_immediate"]     = float(tok[11])
                            ship.info["heightFromSurface"]      = float(tok[12])
                            ship.info["heightFromTerrain"]      = float(tok[13])
                            ship.info["horizontalSrfSpeed"]     = float(tok[14])
                            ship.info["latitude"]               = float(tok[15])
                            ship.info["longitude"]              = float(tok[16])
                            ship.info["obt_velocity"]           = float(tok[17])
                            ship.info["pqsAltitude"]            = float(tok[18])
                            ship.info["rb_velocity"]            = float(tok[19])
                            ship.info["specificAcceleration"]   = float(tok[20])
                            ship.info["srf_velocity"]           = float(tok[21])
                            ship.info["staticPressure"]         = float(tok[22])
                            ship.info["terrainAltitude"]        = float(tok[23])
                            ship.info["verticalSpeed"]          = float(tok[24])
                            ship.info["dynamicPressureAtm"]     = float(tok[25])
                            ship.info["staticPressureAtm"]      = float(tok[26])
                            ship.info["temperature"]            = float(tok[27])
                            
                            
                            ship.ref = toolkit.db.refs[tok[28]]
                            lat = tok[29]
                            lon = tok[30]
                            
                            toolkit.db.active = ship

                            ship.isFlying = False
                            ship.coords = [lat,lon]
                            #print "Landed (%s)"%tok[2]
                            
                        elif tok[1] == "F":
                            pid = tok[2]
                            ship = toolkit.db.vessels[pid]
                            toolkit.db.UT = float(tok[3])
                            toolkit.db.UTt = time.time() 
                            
                            ship.info["missionTime"]            = float(tok[4])
                            ship.info["acceleration"]           = float(tok[5])
                            ship.info["altitude"]               = float(tok[6])
                            ship.info["angularMomentum"]        = float(tok[7])
                            ship.info["angularVelocity"]        = float(tok[8])
                            ship.info["atmDensity"]             = float(tok[9])
                            ship.info["geeForce"]               = float(tok[10])
                            ship.info["geeForce_immediate"]     = float(tok[11])
                            ship.info["heightFromSurface"]      = float(tok[12])
                            ship.info["heightFromTerrain"]      = float(tok[13])
                            ship.info["horizontalSrfSpeed"]     = float(tok[14])
                            ship.info["latitude"]               = float(tok[15])
                            ship.info["longitude"]              = float(tok[16])
                            ship.info["obt_velocity"]           = float(tok[17])
                            ship.info["pqsAltitude"]            = float(tok[18])
                            ship.info["rb_velocity"]            = float(tok[19])
                            ship.info["specificAcceleration"]   = float(tok[20])
                            ship.info["srf_velocity"]           = float(tok[21])
                            ship.info["staticPressure"]         = float(tok[22])
                            ship.info["terrainAltitude"]        = float(tok[23])
                            ship.info["verticalSpeed"]          = float(tok[24])
                            ship.info["dynamicPressureAtm"]     = float(tok[25])
                            ship.info["staticPressureAtm"]      = float(tok[26])
                            ship.info["temperature"]            = float(tok[27])
                            ship.ref = toolkit.db.refs[tok[28]]
                            lat = tok[29]
                            lon = tok[30]
                            
                            toolkit.db.active = ship
                            ship.isFlying = True
                            ship.coords = [lat,lon]
                            #print "Flying (%s)"%tok[2]
                        elif tok[1] == "K":
                            pid = tok[2]
                            ship = toolkit.db.vessels[pid]
                            toolkit.db.UT = float(tok[3])
                            toolkit.db.UTt = time.time() 
                            
                            ship.info["missionTime"]            = float(tok[4])
                            ship.info["acceleration"]           = float(tok[5])
                            ship.info["altitude"]               = float(tok[6])
                            ship.info["angularMomentum"]        = float(tok[7])
                            ship.info["angularVelocity"]        = float(tok[8])
                            ship.info["atmDensity"]             = float(tok[9])
                            ship.info["geeForce"]               = float(tok[10])
                            ship.info["geeForce_immediate"]     = float(tok[11])
                            ship.info["heightFromSurface"]      = float(tok[12])
                            ship.info["heightFromTerrain"]      = float(tok[13])
                            ship.info["horizontalSrfSpeed"]     = float(tok[14])
                            ship.info["latitude"]               = float(tok[15])
                            ship.info["longitude"]              = float(tok[16])
                            ship.info["obt_velocity"]           = float(tok[17])
                            ship.info["pqsAltitude"]            = float(tok[18])
                            ship.info["rb_velocity"]            = float(tok[19])
                            ship.info["specificAcceleration"]   = float(tok[20])
                            ship.info["srf_velocity"]           = float(tok[21])
                            ship.info["staticPressure"]         = float(tok[22])
                            ship.info["terrainAltitude"]        = float(tok[23])
                            ship.info["verticalSpeed"]          = float(tok[24])
                            ship.info["dynamicPressureAtm"]     = float(tok[25])
                            ship.info["staticPressureAtm"]      = float(tok[26])
                            ship.info["temperature"]            = float(tok[27])
                            
                            ship.ref = toolkit.db.refs[tok[28]]
                            epoch = tok[29]
                            a = tok[30]
                            e = tok[31]
                            i = tok[32]
                            W = tok[33]
                            w = tok[34]
                            M = tok[35]
                            
                            toolkit.db.active = ship
                            ship.isFlying = True
                            ship.orbit = toolkit.Orbit(ship,ship.Tref,epoch,a,e,i,W,w,M)
                        else:
                            self.statusmsg.set("Could not parse data")

        try:
            self.toplevel.after_cancel(self.after)
        except:
            pass
        #self.after = self.toplevel.after(1000, self.update)