# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 14:39:54 2013

Phase angle calculator
@author: snaipperi
"""


#import tk
import scipy.optimize as so
from pylab import degrees, arccos, norm, pi, sqrt, floor
import time
#def getSynodicPeriod(from_planet,to_planet):
output = {}
def getTransferAngle(depart_ut, arrive_ut, depart_planet, arrive_planet):
    # Need to project them on the same plane..?
    
    dep = depart_planet.eph(depart_ut)[0]
    dep[2] = 0.0
    #dep_m = norm(dep_v)
    #dep_u = dep_v / dep_m
    dep /= norm(dep)
    arr = arrive_planet.eph(arrive_ut)[0]
    arr[2] = 0.0
    arr /= norm(arr)
    
    #print degrees(arccos(dep.dot(arr)))
    #return o
    return degrees(arccos(dep.dot(arr)))

def get180TransferAngle(depart_ut,depart_planet,arrive_planet):
    synodic = depart_planet.orbit.synodicPeriod(arrive_planet.orbit)
    period = arrive_planet.orbit.period() #+ synodic #Instead we search for the destination whole year..
    #period2 = depart_planet.orbit.period() 
    angle = lambda x: 180 - getTransferAngle(depart_ut,x,depart_planet,arrive_planet)
    boundary =[depart_ut,depart_ut+period]
    #print "Boundary",boundary
    arrive_ut = so.minimize_scalar(angle,method="bounded",bounds=boundary)
    #output[TA.x] = getTransferAngle(depart_ut,TA.x,depart_planet,arrive_planet)
    #print "Optimized at",arrive_ut.x
    #print "Transfer angle",getTransferAngle(depart_ut,arrive_ut.x,depart_planet,arrive_planet)
    
    return arrive_ut.x

def getPeriodDifference(depart_ut,arrive_ut,depart_planet,arrive_planet):
    p1 = abs(depart_ut - arrive_ut)
    
    deph_v = depart_planet.eph(depart_ut)[0]
    aeph_v = arrive_planet.eph(arrive_ut)[0]
    deph = norm(deph_v)
    aeph = norm(aeph_v)
    #print "deph",deph
    #print "aeph",aeph
    #deph_v /= deph
    #aeph_v /= aeph
    #PA = degrees(arccos(deph_v.dot(aeph_v)))
    
    p2 = pi * sqrt((deph+aeph)**3 / (8*depart_planet.ref.mu))
    #print "DEP",depart_ut
    #print "ARR",arrive_ut
    #print "PD",abs(p2-p1)
    #print "P1",p1
    #print "P2",p2
    #print "PA",PA,"degrees"
    #time.sleep(0.1)
    #time.sleep(0.5)
    return abs(p2-p1)

def findNextTransferAngle(depart_ut,depart_planet,arrive_planet):
    
    global output
    output = {}
    
    synodic = depart_planet.orbit.synodicPeriod(arrive_planet.orbit)
    boundary = [depart_ut,depart_ut+synodic]
    
    #get180TransferAngle(depart_ut,depart_planet,arrive_planet)
    
    superfinder = lambda x: getPeriodDifference(x,get180TransferAngle(x,depart_planet,arrive_planet),depart_planet,arrive_planet)
    
    TA = so.minimize_scalar(superfinder,method="bounded",bounds=boundary)
    #print "Start",superfinder(depart_ut)
    #print "End",superfinder(depart_ut+synodic)
    #raw_input("Wait for it")
    #TA = so.brentq(superfinder,depart_ut,depart_ut+synodic)
    #TA = so.brute(superfinder,[boundary])
    #print "We have a solution.."
    #print TA.x
    
    #deph_v = depart_planet.eph(TA.x)[0]
    #aeph_v = arrive_planet.eph(TA.x)[0]
    #deph = norm(deph_v)
    #aeph = norm(aeph_v)
    #print "deph",deph
    #print "aeph",aeph
    #deph_v /= deph
    #aeph_v /= aeph
    #PA = degrees(arccos(deph_v.dot(aeph_v)))    
    #print "PA",PA
    return TA.x
    #pa = lambda x: from_planet.eph(x)[0]
    
def getAnglesBySynodics(synodics,depart_planet,arrive_planet):
    ''' Find all departure dates that have 180 degree transfer angles
    TODO: start searching from non-zero time
    '''    
    
    synodics = int(synodics)
    
    # Get the synodic period
    synodic = depart_planet.orbit.synodicPeriod(arrive_planet.orbit)
    print "The synodic period is",synodic/60.0/60.0/24.0,"days"
    print "Fetching the next",synodics,"synodic periods.."
    print "Which is",synodic/60.0/60.0/24.0*synodics,"days worth of data (",synodic/60.0/60.0/24.0/365.0*synodics," years)"
    
    # Ask for enter so user can confirm that the query matches the requirements
    #raw_input("Hit enter to proceed")
    time.sleep(1)

    # Iterate over all the synodic periods and find a matching 180 degree transfer angle
    departures = []    
    for i in xrange(synodics):
        departures.append(findNextTransferAngle(synodic*i,depart_planet,arrive_planet))
        if i % 10 == 0:
            print "Done",float(i)/synodics*100,"%"
        
    generate_datafile("phaseangle-"+depart_planet.name + "-" + arrive_planet.name + ".txt",departures,depart_planet,arrive_planet)
    return departures
    
    
def getAnglesByYears(years,depart_planet,arrive_planet):
    '''
    Trigger Au, this is the function you're looking for!
    '''
    
    # Get the synodic period
    synodic = depart_planet.orbit.synodicPeriod(arrive_planet.orbit)
    
    # Get the amount of synodic periods during the years
    synodics = (years*365*24*60*60) / synodic
    
    # And proceed
    getAnglesBySynodics(synodics,depart_planet,arrive_planet)
    
def getAngle(t,d,a):
    ''' Just get the angle between two planets at time t '''
    dep = d.eph(t)
    arr = a.eph(t)
    dep = dep[0]
    arr = arr[0]
    dep[2] = 0.0
    arr[2] = 0.0
    dep /= norm(dep)
    arr /= norm(arr)
    return degrees(arccos(dep.dot(arr)))
    

def findAlignment(error,a,b,c):
    ''' 
    find planetary alignments :-)
    error is maximum error in degrees
    a is initial planet
    b is second planet
    c is a list of other planets
    
    choose planets with longest periods as a and b to gain the most efficiency
    also c should be sorted with period length from longest to shortest I suppose
    
    Moho-Kerbin-Eve alignment: 14078521.6823 (Year 1, Day 163)
    Moho-Kerbin-Eve-Duna alignment: 44956459758.5 (Year 1426.. KSP shows negative years)
    Moho-Kerbin-Eve-Duna-Dres-Jool alignment: 2173607115160.0 (Year 68925)
    '''
    synodic = a.orbit.synodicPeriod(b.orbit)
    i = 0
    while True:
        boundary = [i*synodic,(i+1)*synodic]
        f = lambda x: getAngle(x,a,b)
        t = so.minimize_scalar(f,method="bounded",bounds=boundary)
        fa = f(t.x)
        #print "Aligned at",t.x
        #print "Angle",fa
        if fa > error:
            #print "Skipping main, error too large"
            continue
        else:
            fail = False
            for o in c:
                angle = getAngle(t.x,a,o)
                if angle > error:
                    #print "Skipping %s-%s, error too large (%f)"%(a.name,o.name,angle)
                    fail = True
                    break
            if fail:
                i += 1
                if i%100 == 0:
                    print "Current year",t.x / 31536000
                continue
            else:
                print "Planetary alignment found at t",t.x
                break
            
                
        
    
def generate_datafile(name,departures,depart_planet,arrive_planet):
    f = open(name,'w')
    last = departures[-1]
    lasty = int(last / 60.0 / 60.0 / 24.0 / 365.0)
    f.write("// Phase angles for the next %i years for %s - %s Hohmann transfers\n"%(lasty,depart_planet.name,arrive_planet.name))
    f.write("// Calculated using the KSP Mission Control toolkit\n")
    f.write("// Angles are valid for Kerbal Space Program 0.18.4\n\n")
    f.write("UT Departure\tPhase angle\tDate time\n")
    for departure in departures:
        e1 = depart_planet.eph(departure)[0]
        e2 = arrive_planet.eph(departure)[0]
        e1 /= norm(e1)
        e2 /= norm(e2)
        PA = degrees(arccos(e1.dot(e2)))
        
        years = floor(departure/60.0/60.0/24.0/365.0)+1
        days = floor((departure/60.0/60.0/24.0)%365.0)+1
        
        f.write("%f\t%f\tYear %i, Day %i\n"%(departure,PA,years,days))
        
    f.close()
    
    
'''
Use this script to generate stuff for all planets

import toolkit as tk
import phaseangle as pa

planets = [tk.Moho,tk.Eve,tk.Kerbin,tk.Duna,tk.Dres,tk.Jool,tk.Eeloo]
years = 100

for fplanet in planets:
    for tplanet in planets:
        if fplanet == tplanet: continue
        print "Now doing planet",fplanet.name,"to",tplanet.name
        pa.getAnglesByYears(years,fplanet,tplanet)
        
'''