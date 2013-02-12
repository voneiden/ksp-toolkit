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

class Network:
    def __init__(self,ip="192.168.1.13",port=11000):
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        
        self.figure = pl.figure(figsize=(1,2))
        self.axis = self.figure.gca(projection="rectilinear")
        self.buf = ""
        try:        
            self.socket.connect((ip,port))
        except:
            #raise
            self.text = self.axis.text(0.5,0.5,"Network: Failed")
            self.text.set_horizontalalignment("center")
            self.text.set_verticalalignment("center")
            self.figure.canvas.draw()
        else:
            self.socket.setblocking(0)
            self.text = self.axis.text(0.5,0.5,"Network: Connected")
            self.text.set_horizontalalignment("center")
            self.text.set_verticalalignment("center")
            
            self.figure.canvas.manager.window.after(1000, self.update)
            self.figure.canvas.draw()
            
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
                self.text = self.axis.text(0.5,0.5,"Network: Disconnected")
                self.text.set_horizontalalignment("center")
                self.text.set_verticalalignment("center")
                self.figure.canvas.draw()
                return
            else:
                #print "Processing data"
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
                        if tok[1] == "D":
                            pid = tok[2]
                            name = tok[3]
                            #print "Define vessel (%s) %s"%(tok[2],tok[3])
                            toolkit.db.vessels[pid] = toolkit.Ship(name,None)
                            
                        elif tok[1] == "L":
                            pid = tok[2]
                            ref = toolkit.db.refs[tok[3]]
                            lat = tok[4]
                            lon = tok[5]
                            ship = toolkit.db.vessels[pid]
                            
                            toolkit.db.active = ship
                            ship.ref = ref
                            ship.isFlying = False
                            ship.coords = [lat,lon]
                            #print "Landed (%s)"%tok[2]
                            
                        elif tok[1] == "F":
                            pid = tok[2]
                            ref = toolkit.db.refs[tok[3]]
                            lat = tok[4]
                            lon = tok[5]
                            ship = toolkit.db.vessels[pid]
                            
                            toolkit.db.active = ship
                            ship.ref = ref
                            ship.isFlying = True
                            ship.coords = [lat,lon]
                            #print "Flying (%s)"%tok[2]
                        elif tok[1] == "K":
                            pid = tok[2]
                            ref = toolkit.db.refs[tok[3]]
                            epoch = tok[4]
                            a = tok[5]
                            e = tok[6]
                            i = tok[7]
                            W = tok[8]
                            w = tok[9]
                            M = tok[10]
                            ship = toolkit.db.vessels[pid]
                            
                            toolkit.db.active = ship
                            ship.ref = ref
                            ship.isFlying = True
                            ship.orbit = toolkit.Orbit(ship,ref,epoch,a,e,i,W,w,M)
                            #print "Orbiting (%s)"%tok[2]
                
        self.figure.canvas.manager.window.after(1000, self.update)
