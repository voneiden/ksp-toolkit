# -*- coding: utf-8 -*-
"""
Created on Thu Apr 04 21:46:19 2013

@author: snaipperi
"""

''' A second attempt in creating phase angle solver '''
from pylab import degrees, norm, cross, arccos, pi, sqrt, linspace, floor, sin, radians,cos, isnan
import matplotlib.pyplot as plt
import numpy, time

import scipy.optimize as so
STUFF = []

def getAngle(t1,c1,t2,c2):
    ''' 
    Get angle between two celestials at t1 and t2 
    
    Verify if ignoring the k-cordinate makes any sense
    
    timeit 240 microseconds
    '''
    if type(t2) == numpy.ndarray:
            t2 = t2[0]
    elif isnan(t2):
        print "ERROR, t2 is nan!"
        return t2
    
    p1 = c1.eph(t1)[0]
    p1[2] = 0.0
    p1l = norm(p1)
    p1 /= p1l
    
    p2 = c2.eph(t2)[0]
    p2[2] = 0.0
    p2l = norm(p2)
    p2 /= p2l
    
    #if p1l > p2l:
    return p1.dot(p2)
    #else:
    #    return p1.dot(p2)
    
    '''
    if sign:
        
        try:
            c = cross(p1,p2)[2]
        except IndexError:
            print "P1",p1,type(p1)
            print "P2",p2,type(p2)
            global STUFF
            STUFF = [p1,p2,t1,t2]
            raise
        c /= abs(c)
        #angle = degrees(arccos(p1.dot(p2))) * c
        #return angle
        return p1.dot(p2) #* c
    else:
        return degrees(arccos(p1.dot(p2)))
        
        # Could try to skip arccos & degrees to speed up the algorithm
    '''
        
def solveT2(t1,c1,c2,space=11):
    ''' Given t1, solve a valid t2
    base the initial guess on hohmann transfer average
    '''
    #print "SOLVET2 called"
    
    guess = t1 + pi * sqrt((c1.orbit.a+c2.orbit.a)**3 / (8*c1.ref.mu))
    
    ''' Some obsolete sad attempts. Keeping here for reference'''    
    #f1 = lambda t2: 180-getAngle(t1,c1,t2,c2,sign=False)
    #f1 = lambda t2: sin(radians(getAngle(t1,c1,t2,c2,sign=True)))
    # Don't use sin!, instead cos and divide the angle by two
    #f1 = lambda t2: cos(radians(getAngle(t1,c1,t2,c2,sign=True)) / 2.0)
    
    
    f1 = lambda t2: sqrt(getAngle(t1,c1,t2,c2) + 1) / 1.4142135623730951
    # Unoptimized version: cos(arccos(x)/2.0)
    # Optimized version  : sqrt(x+1) / sqrt(2) (=1.4142135623730951)
    
    try:
        '''
        plt.figure()
        plotspace = linspace(t1,guess+(guess-t1),50)
        values = []

        for s in plotspace:
            print "Trying space",s
            values.append(f1(s))

        plt.plot(plotspace,values)
        '''

        
        # Setup initial conditions
        value_guess = f1(guess)
        
        """
        print "Check"
        print t1,type(t1)
        print guess,type(guess)
        print space,type(space)
        """
        searchspace = linspace(t1,guess+(guess-t1),space)
        
        i = (space-1) / 2 # Default 25
        
        value_up = f1(searchspace[i+1])
        value_down = f1(searchspace[i-1])
        
        up_decreasing = value_up < value_guess
        down_decreasing = value_down < value_guess
        
        # This situation happens if the initial guess leads to a solution
        # that's either very close to zero or very close to 1
        # Neither direction looks like it has the closest correct solution in..
        if (up_decreasing and down_decreasing) or (not up_decreasing and not down_decreasing):
            '''           
            print "Solver is stalled!!"
            print "Guess value:",value_guess
            print "Up value:   ",value_up
            print "Down value: ",value_down
            '''
            # If it's very close to zero, we have it easy. We can just bound
            # the space around the initial guess and get the value out
            if value_guess < 0.5:
                t2 = so.minimize_scalar(f1,method="Bounded",bounds=[searchspace[i-1],searchspace[i+1]]).x
            
            # On the other hand, if the value is close to one, 
            # we need either more accuracy or just choose an arbitary
            # direction (I believe this causes a discontinuity in the f1 function)            
            else:
                #print "INCREASING SPACE"
                t2 = solveT2(t1,c1,c2,space=(space-1)*10+1)
            
            
        else:
            # We tried increasing the space (51.. 501.. 5001..)
            # But 501 didn't lead to any solution either.
            # Lets get arbitary, and just.. go up.
            if space > 1000:
                #print "Space is over 1000!"
                if value_guess > 0.5:
                    # Lets just search up
                    #print "Looks like 0.999999... Searching up for the sake of humanity"
                    
                    #This could be optimized a bit
                    space = 11
                    i = 10
                    searchspace = linspace(t1,guess+(guess-t1),space)
                    up_decreasing = True
                else:
                    print "Error, unable to continue. Guess value is close to zero but does not solve."
                    raise RuntimeError
                       
                        
            search_up = None
            
            if up_decreasing:
                search_up = True
            else:
                search_up = False
            
            '''
            print "Guess value:",value_guess
            print "Up value:   ",value_up
            print "Down value: ",value_down
            print "Search up? :",search_up
            '''
            
            last_value1 = [value_guess,searchspace[i]]
            last_value2 = [value_guess,searchspace[i]]    
            
            while True:
                if search_up:
                    i += 1
                else:
                    i -= 1
                    if i < 0:
                        #print "Ran out of bounds while going down.. enlarging search area.."
                        #diff = searchspace[1]-searchspace[0]
                        #searchspace = numpy.insert(searchspace,0,searchspace[0]-diff)
                        #i = 0
                        
                        # Optimization attempt
                        boundary = searchspace[-1] - searchspace[0]
                        #diff = searchspace[1] - searchspace[0]
                        
                        searchspace = linspace(searchspace[1]-boundary,searchspace[1],space)
                        i = len(searchspace) - 2
                        continue
                        
                        
                        
    
                try:
                    value = f1(searchspace[i])
                except IndexError:
                    #print "Ran out of bounds while going up.. need to enlarge the search area"
                    
                    if search_up:
                        #diff = searchspace[i-1]-searchspace[i-2]
                        #searchspace = numpy.append(searchspace, searchspace[i-1]+diff)
                        #i -= 1
                        """
                        print "Diff is",diff
                        print "New max is",searchspace[i-1]+diff
                        """
                        
                        boundary = searchspace[-1] - searchspace[0]
                        #diff = searchspace[1] - searchspace[0]
                        """
                        print "Current space",searchspace[0],searchspace[-1]
                        print "The boundary is",boundary
                        print "Going for a searchspace",searchspace[-2],searchspace[-2]+boundary
                        """
                        searchspace = linspace(searchspace[-2],searchspace[-2]+boundary,space)
                        i = 1
                        
                        continue
                    else:
                        print "FAIL"
                        raise RuntimeError
                   
                #print "CHECK"
                if value > last_value1[0]:
                    if search_up:
                        boundary = [last_value2[1],searchspace[i]]
                    else:
                        boundary = [searchspace[i],last_value2[1]]
                    """
                    print "value",value
                    print "lastvalue1",last_value1
                    print "lastvalue2",last_value2
                    """
                    break
                else:
                    last_value2 = [last_value1[0],last_value1[1]]
                    last_value1 = [value,searchspace[i]]
                    
            #print "Solution found in boundary:",boundary
            #try:
            t2 = so.minimize_scalar(f1,method="Bounded",bounds=boundary).x
            #except ValueError:
            #    print "Attempted boundary is",boundary
            #    raise
            


       # print "Root found at",t2
        #if f1(t2) > 1:
        #    print "The root has an angle of",f1(t2)
        #    time.sleep(1)
        return t2
        
        
        
        
        
    except:
        print "Something blew up. Good luck debugging."
        #raise
        '''
        print "Newton method failed, why? Guess:",guess
        plt.figure()
        space = linspace(t1,guess+(guess-t1),50)
        values = []

        for s in space:
            values.append(f1(s))

        plt.plot(space,values)
        '''
        raise
        
        
        
        """
        raise
        print "Try newton again"
        f2 = lambda t2: getAngle(t1,c1,t2,c2,sign=True)
        try:
            t2 = so.newton(f2,guess,tol=0.1,maxiter=500)
        except:
            plt.figure()
            space = linspace(t1,guess+(guess-t1),50)
            values = []

            for s in space:
                print "Trying space",s
                values.append(f1(s))

            plt.plot(space,values)
            
            plt.figure()
            space = linspace(t1,guess+(guess-t1),50)
            values = []

            for s in space:
                values.append(f2(s))

            plt.plot(space,values)
            print "RAISING HELL"
            raise
            
        return t2
        """
        

def solveToF(t1,c1,t2,c2,sign=False):
    if t2 == False:
        return False
    tf1 = t2-t1
    
    c1p = c1.eph(t1)[0]
    c2p = c2.eph(t2)[0]
    
    #c1p[2] = 0.0
    #c2p[2] = 0.0
    
    c1p = norm(c1p)
    c2p = norm(c2p)
    
    tf2 = pi * sqrt((c1p+c2p)**3 / (8*c1.ref.mu))
    #if abs(tf1-tf2) > 3000000:
    #    print "OKAY WIERD SITUATION"
    #    print
    #    print
    #    print
    #    print "TF1",tf1
    #    print "TF2",tf2
    if sign:
        return tf1-tf2
    else:
        return abs(tf1-tf2)

def solveT1(boundary,c1,c2,sign=False):
    
    f = lambda t1: solveToF(t1, c1, solveT2(t1,c1,c2), c2,sign)
    
    #try:        
    #    t1 = so.newton(f,guess)
    #except:
    #    print "T1 Newton method failed, attempting scalar"  
    #    f1 = lambda t2: getAngle(t1,c1,t2,c2)
    
    space = linspace(boundary[0],boundary[1])
    # Find the lowest points
    min1 = [999999999,None,0]
    #min2 = [999999999,None]
    
    for i,point in enumerate(space):
        x = f(point)
        if x < min1[0]:
            #min2 = min1
            min1 = [x,point,i]
        #elif x < min2[0]:
        #    min2 = [x,point]
            
    
    
    try:
        if min1[2] == 0:
            raise IndexError
        left = space[min1[2]-1]
    except IndexError:
        left = boundary[0]
    try:
        right = space[min1[2]+1]
    except IndexError:
        right = boundary[1]
    miniboundary = [left,right]
    
    t1 = so.minimize_scalar(f,method="Bounded",bounds=boundary).x
    
    t1 = so.minimize_scalar(f,method="Bounded",bounds=miniboundary).x
    
    t2 = solveT2(t1,c1,c2)
    check = solveToF(t1,c1,t2,c2)
    
    
    
    """
    space= linspace(boundary[0],boundary[1],100)
    
    values = []
    for s in space:
        values.append(f(s))
        
    plt.plot(space,values)
    """
    if check > 100:
        print "There is some error in the check value, investigate?"
        print "Check error:",check
        print "Solution v1:",t1
        print "Solution v2:",t1
        if check > 10000:
            raise RuntimeError
        else: 
            return (t1,check)
    
    """
    print "Boundary    :",boundary
    print "Miniboundary:",miniboundary
    print "Solution v1:",t1
    print "Solution v2:",t1
    print "Check",check
    print "T1 solved at",t1
    """
    return (t1,0)
    
    
    
    
def test(dc,ac,sign=False):
    synodic = dc.orbit.synodicPeriod(ac.orbit)
    plt.figure()
    solveT1([0,synodic],dc,ac,sign)
    
    

def getSynodic(synodicN,depart_planet,arrive_planet,sign=False):
    synodic = depart_planet.orbit.synodicPeriod(arrive_planet.orbit)
    #plt.figure()
    mod = getSynodicModifier(depart_planet,arrive_planet)
    for i in xrange(1):
        t = (synodicN+i)*synodic+mod
        solveT1([t,t+synodic],depart_planet,arrive_planet,sign)
        
def getSynodicModifier(c1,c2):
    f = lambda t1: - solveToF(t1, c1, solveT2(t1,c1,c2), c2)
    synodic = c1.orbit.synodicPeriod(c2.orbit)
    print "Bounds",[-synodic,synodic]
    mod = so.minimize_scalar(f,method="Bounded",bounds=[-synodic,synodic])
    print "Modifier:",mod.x
    return mod.x
    
    
def getBySynodics(synodics,depart_planet,arrive_planet):
    synodic = depart_planet.orbit.synodicPeriod(arrive_planet.orbit)
    mod = getSynodicModifier(depart_planet,arrive_planet)
    
    print "%s - %s transfer"%(depart_planet.name,arrive_planet.name)
    print "The synodic period is",synodic/60.0/60.0/24.0,"days"
    print "Fetching the next",synodics,"synodic periods.."
    print "Which is",synodic/60.0/60.0/24.0*synodics,"days worth of data (",synodic/60.0/60.0/24.0/365.0*synodics," years)"
    
    if synodics == 0:
        print "Skipping since no synodic periods"
        return []
    # Ask for enter so user can confirm that the query matches the requirements
    #raw_input("Hit enter to proceed")
    time.sleep(1)

    # Iterate over all the synodic periods and find a matching 180 degree transfer angle
    departures = []    
    #plt.figure()
    for i in xrange(synodics):
        syn = synodic*i+mod
        print "Synodic",i
        depart_UT = solveT1([syn,syn+synodic],depart_planet,arrive_planet,False)
        departures.append(depart_UT)
        if i % 10 == 0:
            print "Done",float(i)/synodics*100,"%"
        #if TA == False:
        #    print "Error encountered, aborting"
        #    errors += 1
        #else:
        #    departures.append(TA)
    print "Completed. Generating datafile"
    generate_datafile("phaseangle-"+depart_planet.name + "-" + arrive_planet.name + ".txt",departures,depart_planet,arrive_planet)
    print "Datafile generated"    
    return departures
    
    
def getByYears(years,depart_planet,arrive_planet):
    '''
    Trigger Au, this is the function you're looking for!
    '''
    
    # Get the synodic period
    synodic = depart_planet.orbit.synodicPeriod(arrive_planet.orbit)
    
    # Get the amount of synodic periods during the years
    synodics = int((years*365*24*60*60) / synodic)
    
    # And proceed
    print "Calling getBySynodics"
    getBySynodics(synodics,depart_planet,arrive_planet)
    print "GetBySynodics call completed"
    
def getByYearsPlanet(years,depart_planet,planets):
    for planet in planets:
        if planet == depart_planet:
            continue
        getByYears(years,depart_planet,planet)
    
    
def generateEverything(years,planets):
    for dep in planets:
        for arr in planets:
            if dep == arr:
                continue
            else:
                getByYears(years,dep,arr)
                
def generate_datafile(name,departures,depart_planet,arrive_planet,errors=0):
    f = open(name,'w')
    last = departures[-1]
    lasty = int(last[0] / 60.0 / 60.0 / 24.0 / 365.0)
    inacc = 0.0
    for departure in departures:
        if departure[1]:
            inacc += 1.0
            
    f.write("// Phase angles for the next %i years for %s - %s Hohmann transfers\n"%(lasty,depart_planet.name,arrive_planet.name))
    f.write("// Calculated using the KSP Mission Control toolkit\n")
    f.write("// Angles are valid for Kerbal Space Program 0.19\n")
    f.write("// Total windows: %i\n"%len(departures))
    f.write("// Inaccuracies during calculation: %i (%i%%)\n\n"%(inacc,inacc / len(departures) * 100))
    f.write("UT Departure\tPhase angle\tDate time\tAccurate (Error in seconds)\n")
    for departure in departures:
        accurate = departure[1]
        if accurate == 0:
            accurate = "Yes"
        else:
            accurate = str(accurate)
            
        departure = departure[0]
        
        e1 = depart_planet.eph(departure)[0]
        e2 = arrive_planet.eph(departure)[0]
        e1 /= norm(e1)
        e2 /= norm(e2)
        PA = degrees(arccos(e1.dot(e2)))
        
        years = floor(departure/60.0/60.0/24.0/365.0)+1
        days = floor((departure/60.0/60.0/24.0)%365.0)+1
        
        f.write("%f\t%f\tYear %i, Day %i\t%s\n"%(departure,PA,years,days,accurate))
        
    f.close()