# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

'''QUESTIONS:'''
# 1. uu is selected correctly. But is dd? I guess so, but not sure...
#    EV TODO = try it differently. Does is work now/still?

'''TODO'''
# 1.

'''PROGRAM'''
class CurveDesigner(object):
    """Class to design 2D curves using cubic splines and control points, 
    as well as by fitting cubic splines through given data points."""

    def __init__(self):
        pass

    def deBoor(self, dd, ui, u):
        """"Calculates new points s(u) on a curve using the De Boor algorithm.
        dd: [d_(I-2), ..., d_(I+1)] : control points.
        ui: [u_(I-2), ..., u_(I+3)] : node points.
        u: parameter for generating a new point."""
       
        uu2, uu1, u0, u1, u2, u3 = ui
        
        a21 = (u1-u)/(u1-uu2)
        a31 = (u2-u)/(u2-uu1)
        a41 = (u3-u)/(u3-u0)
        d21 = a21*dd[0] + (1-a21)*dd[1]
        d31 = a31*dd[1] + (1-a31)*dd[2]
        d41 = a41*dd[2] + (1-a41)*dd[3]
        hfsejkfhsejk = 1
        a32 = (u1-u)/(u1-uu1)
        a42 = (u2-u)/(u2-u0)
        d32 = a32*d21 + (1-a32)*d31
        d42 = a42*d31 + (1-a42)*d41
        
        a43 = (u1-u)/(u1-u0)
        d43 = a43*d32 + (1-a43)*d42
        
        return d43

#Define the grid and the control points.
d_vector = np.array([[1, 4], [0.5, 6], [5, 4], [3, 12], [11, 14], 
      [8, 4], [12, 3], [11, 9], [15, 10], [17, 8]]) #initial control points (10)
u_vector = np.array([0,0,0,0,1/8,2/8,3/8,4/8,5/8,6/8,1,1,1,1]) #knot vector (14)

#Find the "hot interval".
u = 0.99
i = int(u_vector.searchsorted([u])) #6 for u=0.3.

#Select the six relevant u's and the four relevant d's.
uu = u_vector[i-3:i+3] #
dd = d_vector[i-4:i]                                                            # Correct?

#Blossom recursion.
cd = CurveDesigner()
new_u = cd.deBoor(dd, uu, u)
# but we need many new points... hmmm... and then to plot them.






#Plotting control points.
d0, d1 = zip(*d_vector)
plt.plot(d0, d1, 'bo')
# TODO = (1) make these points red instead. (2) draw lines in between them.
    







