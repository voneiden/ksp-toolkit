ksp-toolkit
===========

Dependencies: Python 2.7, matplotlib and pylab (scipy+numpy)
Recommended: ipython (ipythonqt)

Sample commands:

import toolkit as tk
o3d = tk.Orbit3D(tk.Kerbol)
o3d.track(tk.Kerbin)
o3d.track(tk.Duna)
o3d.track(tk.Eve)
o3d.track(tk.Moho)
o3d.track(tk.Jool)

o2d = tk.Orbit2D() # 2d tracks ignore rotation completely 
o2d.track(tk.Mun) # So it only shows the orbit shape from 2d point of view

ts = tk.testShip() # Creates a debug test ship on Kerbin orbit
o2d.track(ts)