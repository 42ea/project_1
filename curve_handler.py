# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

import unittest
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
        if interpolation_points is None:
            pass
        else:
            self.interpolation_points = interpolation_points
            self.xi = (u_vector[:-2]+u_vector[1:-1]+u_vector[2:])/3
            
    
    def __call__(self, d_vector = None, u_vector = None): #Method for using a created CurveDesigner-instance
        
        if d_vector is not None: 
            self.d_vector = d_vector

        if u_vector is not None: 
            self.u_vector = u_vector
            
    def generateSpline(self, n, mode):
        ''' This function takes in control points (d_vector) and node points 
        (u_vector) and returns the cubic spline for those points using the 
        deBoor method for cubic splines
        
        n determines the "resolution" of the spline
        
        '''
        if mode=='interpol banded':
            L2 = len(self.u_vector)-2
            Ni = (L2)*[None]
            ab = np.zeros((4,L2))
            for j in range(L2):
                Ni[j] = self.basis_func(j)
                if j==0:
                    ab[0][j] = 0
                    ab[1][j] = 0
                    ab[2][j] = Ni[j](self.xi[j])
                    ab[3][j] = Ni[j](self.xi[j+1])
                if j==1:
                    ab[0][j] = 0
                    ab[1][j] = Ni[j](self.xi[j-1])
                    ab[2][j] = Ni[j](self.xi[j])
                    ab[3][j] = Ni[j](self.xi[j+1])
                if (j>1) and (j<L2-1):
                    ab[0][j] = Ni[j](self.xi[j-2])
                    ab[1][j] = Ni[j](self.xi[j-1])
                    ab[2][j] = Ni[j](self.xi[j])
                    ab[3][j] = Ni[j](self.xi[j+1])
                if j==L2-1:
                    ab[0][j] = Ni[j](self.xi[j-2])
                    ab[1][j] = Ni[j](self.xi[j-1])
                    ab[2][j] = Ni[j](self.xi[j])
                    ab[3][j] = 0
        if mode == 'interpol not banded':
            L2 = len(self.u_vector)-2
            Ni = (L2)*[None]
            A = np.zeros((L2,L2))
            ab = np.zeros((4,L2))
            for j in range(L2):
                Ni[j] = self.basis_func(j)
                for i in range(L2):
                    A[i][j] = Ni[j](self.xi[i])
            d_x = np.linalg.solve(A, self.interpolation_points[:,0])
            d_y = np.linalg.solve(A, self.interpolation_points[:,1])
            CONTROL_POINTS = np.column_stack((d_x, d_y))
            self.d_vector=CONTROL_POINTS

        elif mode=='control points':
            self.u = np.linspace(min(self.u_vector)+0.001,max(self.u_vector)-0.001,n)  # generate n- long vector of u-values to generate spline
            spline = np.empty([2,n])
            
            #This for-loop can probably be replaced with vector operations 
            for j in range(0,n):
                i = int(self.u_vector.searchsorted([self.u[j]]))   # finds the "hot interval"
                uu = self.u_vector[i-3:i+3] #extracts the relevant 
                dd = self.d_vector[i-3:i+1]   
                S_u = self.deBoor(dd,uu,self.u[j]) # generate the S value for our current u. 
                spline[:,j] = S_u  # add the the point to the spline vector
            return spline
        pass
    
    
    
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
            #Prevent out of bounds and check if divide by zero.
            if j==0:
                if (self.u_vector[j+k]-self.u_vector[j]):
                    return (self.u_vector[j+k]-u)/(self.u_vector[j+k]-self.u_vector[j])*N(u,j+1,k-1)
                else:
                    return 0
            
            if (j+1)==(len(self.u_vector)-2):
                if (self.u_vector[j+k-1]-self.u_vector[j-1]):
                    return (u-self.u_vector[j-1])/(self.u_vector[j+k-1]-self.u_vector[j-1])*N(u,j,k-1)    
                else:
                    return 0

            #Not out of bounds, check if divide by zero
            if ((self.u_vector[j+k-1]-self.u_vector[j-1]) and (self.u_vector[j+k]-self.u_vector[j])):
                return (u-self.u_vector[j-1])/(self.u_vector[j+k-1]-self.u_vector[j-1])*N(u,j,k-1) \
                        +(self.u_vector[j+k]-u)/(self.u_vector[j+k]-self.u_vector[j])*N(u,j+1,k-1)
            
            if ((self.u_vector[j+k-1]-self.u_vector[j-1]) and not (self.u_vector[j+k]-self.u_vector[j])):
                return (u-self.u_vector[j-1])/(self.u_vector[j+k-1]-self.u_vector[j-1])*N(u,j,k-1)
            
            if ((self.u_vector[j+k]-self.u_vector[j]) and not (self.u_vector[j+k-1]-self.u_vector[j-1])):
                return (self.u_vector[j+k]-u)/(self.u_vector[j+k]-self.u_vector[j])*N(u,j+1,k-1)
        
        def evaluate_N(u):
            return N(u,j,3) 
        return evaluate_N
    
    def plot(self, spline, d_vector, control = False):
        s1 = spline[0,:]        # generate the x-coordinates for the spline
        s2 = spline[1,:]        # generate the y-coordniates for the spline

        d0, d1 = zip(*d_vector) #Separates the x- and y-values for the control points
        
        if control:
            #Plot control points with line inbetween
            plt.plot(d0, d1, color = 'r', linewidth = 0.3)
            plt.plot(d0, d1, 'bo', color = 'r')
        plt.plot(s1,s2)         # plotting the spline 
        plt.show()
        
    def splineFromBasisFunc(self, n):
        #This function evaluates S(u) for each u using the basis functions 
        #given from basis_func
        
        #Create empty list for storing L basis functions
        Ni = (len(self.u_vector)-2)*[None]
        
        #Create the basis functions for each value j (each possible hot interval) and store them in a list
        for j in range(len(self.u_vector)-2):
            Ni[j] = self.basis_func(j)
            
        #Use the control points and the basis functions to create s(u)
        
        self.u = np.linspace(min(self.u_vector)+0.001,max(self.u_vector)-0.001,n)  #Generate n-vector of u-values to generate spline
        Spline = np.empty([2,n]) #Create empty vector for filling with points on the spline
        
        #For each u in the u-vector, compute the spline and insert into spline vector
        for j in range(0,n):
            i = int(self.u_vector.searchsorted([self.u[j]])) #Find the "hot interval"
            
            dd = self.d_vector[i-3:i+1] #Extract the useful values of d_i from i = I - 2 to i = I + 1
            Nfuncs = Ni[i-3:i+1] #Extract the useful values of N^3_i from i = I - 2 to i = I + 1
            
            #Generate the spline for our current u by summing the products of each d_i with its corresponding N^3_i
            S_u = np.empty([])
            for p in range(len(dd)):
                Nval = Nfuncs[p](self.u[j])
                S_u = S_u + dd[p, :]*Nval
                
            Spline[:,j] = S_u  # add the the point to the spline vector   
        return Spline
            
        