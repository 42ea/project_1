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

    def __init__(self, d_vector = None, u_vector = None, interpolation_points = None): #Constructor method for a CurveDesigner-object
        
        if d_vector is None:
            self.d_vector = np.array([[1, 4], [0.5, 6], [5, 4], [3, 12], [11, 14], 
                                            [8, 4], [12, 3], [11, 9], [15, 10], [17, 8]]) #Default control point-vector
        else:
            self.d_vector = d_vector
            
        if u_vector is None: 
            self.u_vector = np.array([0, 0, 0, 0, 1/8, 2/8, 3/8, 4/8, 5/8, 6/8, 1, 1, 1, 1]) #Default knot vector
        else:
            self.u_vector = u_vector
            
    
    def __call__(self, d_vector = None, u_vector = None): #Method for using a created CurveDesigner-instance
        
        if d_vector is not None: 
            self.d_vector = d_vector

        if u_vector is not None: 
            self.u_vector = u_vector
            
    def generateSpline(self, d_vector ,u_vector, n):
        ''' This funktion takes in control points (d_vector) and node points 
        (u_vector) and returns the cubic spline for those points using the 
        deBoor method for cubic splines
        
        n determines the "resolution" of the spline
        
        '''
        self.u_vector = u_vector     #
        self.d_vector = d_vector
        self.u = np.linspace(min(self.u_vector)+0.001,max(self.u_vector)-0.001,n)  # generate n- long vector of u-values to generate spline
        Spline = np.empty([2,n])
        
        #This for-loop can probably be replaced with vector operations 
        for j in range(0,n):
            i = int(self.u_vector.searchsorted([self.u[j]]))   # finds the "hot interval"
            uu = self.u_vector[i-3:i+3] #extracts the relevant 
            dd = self.d_vector[i-4:i]   
            S_u = self.deBoor(dd,uu,self.u[j]) # generate the S value for our current u. 
            Spline[:,j] = S_u  # add the the point to the spline vector
        return Spline
    
    
    
    def deBoor(self, dd, ui, u):
        """"Calculates new points s(u) on a curve using the De Boor algorithm.
        dd: [d_(I-2), ..., d_(I+1)] : control points for our hot interval
        uu: [u_(I-2), ..., u_(I+3)] : node points for our 
        u: parameter for generating a new point."""
        uu2, uu1, u0, u1, u2, u3 = ui
        
        a21 = (u1-u)/(u1-uu2)
        a31 = (u2-u)/(u2-uu1)
        a41 = (u3-u)/(u3-u0)
        d21 = a21*dd[0] + (1-a21)*dd[1]
        d31 = a31*dd[1] + (1-a31)*dd[2]
        d41 = a41*dd[2] + (1-a41)*dd[3]
        
        a32 = (u1-u)/(u1-uu1)
        a42 = (u2-u)/(u2-u0)
        d32 = a32*d21 + (1-a32)*d31
        d42 = a42*d31 + (1-a42)*d41
        
        a43 = (u1-u)/(u1-u0)
        d43 = a43*d32 + (1-a43)*d42
        
        return d43
        
    def basis_func(self, j):
        def N(u,j,k):
            if(k==0):
                if(self.u_vector[j]==self.u_vector[j-1]):
                    return 0
                if((self.u_vector[j-1]<=u) and (u<self.u_vector[j])):
                    return 1
                return 0

            return (u-u[j-1])/(u[j+k-1]-u[j-1])*N(u,j,k-1) \
                    +(u[j+k]-u)/(u[j+k]-u[j])*N(u,j+1,k-1)
        def evaluate_N(u):
            return N(u,j,3)
        return evaluate_N
    
    def plot(self, Spline, d_vector, control = False):
        s1 = spline[0,:]        # generate the x-coordinates for the spline
        s2 = spline[1,:]        # generate the y-coordniates for the spline

        d0, d1 = zip(*d_vector) #Separates the x- and y-values for the control points
        
        if control:
            #Plot control points with line inbetween
            plt.plot(d0, d1, color = 'r', linewidth = 0.3)
            plt.plot(d0, d1, 'bo', color = 'r')
            
        plt.plot(s1,s2)         # plotting the spline 
            
        



#Blossom recursion.
cd = CurveDesigner()

#Find the "hot interval".
u = 0.01 # om man börjar u- på noll så kommer uu bli tomt set så frågan är hur vi ska välja första u
i = int(cd.u_vector.searchsorted([u])) #6 for u=0.3.

#Select the six relevant u's and the four relevant d's.
uu = cd.u_vector[i-3:i+3] #
dd = cd.d_vector[i-4:i]                                                            # Correct?

new_u = cd.deBoor(dd, uu, u)
# but we need many new points... hmmm... and then to plot them.
spline = cd.generateSpline(cd.d_vector, cd.u_vector,50)
   
#Plot s(u) and control points
cd.plot(spline, cd.d_vector, control = True, )

# TODO = (1) make these points red instead. (2) draw lines in between them.


