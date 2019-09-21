# -*- coding: utf-8 -*-
import numpy as np
import scipy.linalg as sl
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

class CurveDesigner(object):
    """Class to design curves using cubic splines and control points, 
    as well as by fitting cubic splines through given data points."""

    def __init__(self, d_vector=None, u_vector=None, data_points=None): #Constructor method for a CurveDesigner-object
        
        if d_vector is None:
            self.d_vector = self.load_default_data('d_vector') #Default control point-vector
        else:
            self.d_vector = d_vector
            
        if u_vector is None: 
            self.u_vector = self.load_default_data('u_vector') #Default knot vector
            print(self.u_vector)
        else:
            self.u_vector = u_vector
        if data_points is None:
            pass
        else:
            self.data_points = data_points
            self.xi = (self.u_vector[:-2]+self.u_vector[1:-1]+self.u_vector[2:])/3
            self.xi[-1]=self.xi[-1]-0.001*self.xi[-1]
            
    
    def __call__(self, d_vector = None, u_vector = None,data_points = None): #Method for using a created CurveDesigner-instance
        
        if d_vector is not None: 
            self.d_vector = d_vector

        if u_vector is not None: 
            self.u_vector = u_vector
        if data_points is not None:
            self.data_points = data_points
            self.xi = (u_vector[:-2]+u_vector[1:-1]+u_vector[2:])/3
            self.xi[-1]=self.xi[-1]-0.001*self.xi[-1]
            
            
    def git_test(self):
        print("Please God, let me succeed in pushing to git. I want it soo much.")
        
        
        
    def generate_spline(self, n, interpolate=False):
        ''' This function takes in control points (d_vector) and node points 
        (u_vector) and returns the cubic spline for those points using the 
        deBoor method for cubic splines.
        The parameter n determines the "resolution" of the spline.
        
        If interpolate is equal to true, the method will calculate the control
        points through an interpolation method and then calculate the spline s(u)
        with the newly acquired control points.
        '''
        if interpolate:
            L2 = len(self.u_vector)-2
            Ni = (L2)*[None]
            ab = np.zeros((5,L2))
            for j in range(L2):
                Ni[j] = self.basis_func(j)
                if j==0:
                    ab[0][j] = 0
                    ab[1][j] = 0
                    ab[2][j] = Ni[j](self.xi[j])
                    ab[3][j] = Ni[j](self.xi[j+1])
                    ab[4][j] = Ni[j](self.xi[j+2])
                if j==1:
                    ab[0][j] = 0
                    ab[1][j] = Ni[j](self.xi[j-1])
                    ab[2][j] = Ni[j](self.xi[j])
                    ab[3][j] = Ni[j](self.xi[j+1])
                    ab[4][j] = Ni[j](self.xi[j+2])
                if (j>1) and (j<L2-2):
                    ab[0][j] = Ni[j](self.xi[j-2])
                    ab[1][j] = Ni[j](self.xi[j-1])
                    ab[2][j] = Ni[j](self.xi[j])
                    ab[3][j] = Ni[j](self.xi[j+1])
                    ab[4][j] = Ni[j](self.xi[j+2])
                if j==L2-2:
                    ab[0][j] = Ni[j](self.xi[j-2])
                    ab[1][j] = Ni[j](self.xi[j-1])
                    ab[2][j] = Ni[j](self.xi[j])
                    ab[3][j] = Ni[j](self.xi[j+1])
                    ab[4][j] = 0
                if j==L2-1:
                    ab[0][j] = Ni[j](self.xi[j-2])
                    ab[1][j] = Ni[j](self.xi[j-1])
                    ab[2][j] = Ni[j](self.xi[j])
                    ab[3][j] = 0
                    ab[4][j] = 0
            d_x = solve_banded((2,2),ab, self.data_points[:,0])
            d_y = solve_banded((2,2),ab, self.data_points[:,1])
            self.d_vector=np.column_stack((d_x, d_y))
            
            print('Test norm(ab) == 1')
            print(sl.norm(ab,np.inf))
            print('Test Ni[0](0.001) == 1')
            print(Ni[0](0.001))

        # Create points u, used for plotting the spline S(u)
        self.u = np.linspace(min(self.u_vector)+0.001,max(self.u_vector)-0.001,n) 
        spline = np.empty([np.shape(self.d_vector)[1],n])

        for j in range(0,n):
            i = int(self.u_vector.searchsorted([self.u[j]]))   # Finds the hot interval
            uu = self.u_vector[i-3:i+3] # Extract the relevant interval of knot points  
            dd = self.d_vector[i-3:i+1] # Extract the relevant interval of control points  
            S_u = self.deBoor(dd,uu,self.u[j]) # Generate the S value for our current u. 
            spline[:,j] = S_u  # Add the the point to the spline vector
        return spline
    
    
    
    def deBoor(self, dd, ui, u, get_blossoms = False):
        """"Calculates new points s(u) on a curve using the De Boor algorithm.
        dd: [d_(I-2), ..., d_(I+1)] : control points for our hot interval
        uu: [u_(I-2), ..., u_(I+3)] : node points for our hot interval
        ui: node point array
        u: parameter for generating a new point.
        blossom: parameter for choosing whether
        or not to plot blossom points

        blossom recursion:
            dd[0]
            dd[1]   d21
            dd[2]   d31    d32
            dd[3]   d41    d42    d43

        """
        uu2, uu1, u0, u1, u2, u3 = ui
        
        # Calculating the relevant alpha values and blossoms for the first iteration.
        a21 = (u1-u)/(u1-uu2)
        a31 = (u2-u)/(u2-uu1)
        a41 = (u3-u)/(u3-u0)
        d21 = a21*dd[0] + (1-a21)*dd[1]
        d31 = a31*dd[1] + (1-a31)*dd[2]
        d41 = a41*dd[2] + (1-a41)*dd[3]
        
        # Second iteration
        a32 = (u1-u)/(u1-uu1)
        a42 = (u2-u)/(u2-u0)
        d32 = a32*d21 + (1-a32)*d31
        d42 = a42*d31 + (1-a42)*d41
        
        # Third and final iteration
        a43 = (u1-u)/(u1-u0)
        d43 = a43*d32 + (1-a43)*d42
        
        #If blossom is set to True, return the blossoms, from the first and second iteration, and control point for u_blossom
        #instead of the spline point
        if get_blossoms:
            blossoms = d21, d32 #Defines the second iteration blossom points d[u, u, u_I] 
            control_point = dd[0] #Defines the control point
            return blossoms, control_point, d43
        
        return d43
        
    def basis_func(self, j):
        """
            By specifying a grid point j, one will obtain a function
            handle that will be able to evaluate basis function j at a specified u-value.
            in: j
            out: N_j(u)
        """
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
    
    def plot(self, spline, d_vector, control=False, interpolate=False, blossom=False, blossoms=None, control_point=None, d43=None):
        s1 = spline[0,:]        # generate the x-coordinates for the spline
        s2 = spline[1,:]        # generate the y-coordniates for the spline
        plt.plot(s1,s2)         # plotting the spline 
        
        if control:
            d0, d1 = zip(*d_vector) #Separates the x- and y-values for the control points
            plt.plot(d0, d1, color = 'r', linewidth = 0.3) #Plot control points with line inbetween
            plt.plot(d0, d1, 'bo', color = 'r')
        
        if blossom: # Plot ALL blossom points for creating a certain S(u), instead of the two required blossom points in add-on 1
            b = np.array(blossoms)
            plt.plot(d43[0], d43[1], 'bo', color = 'b') #Plot the given spline point in blue.
            plt.plot(b[0,0], b[0,1], 'X', color = 'g', markersize = 8) #Plot the first iteration blossom point d21 = d[u, u_{I-1}, u_I] as a green cross.
            plt.plot(b[1,0], b[1,1], 'X', color = 'y', markersize = 8) #Plot the second iteration blossom point d32 = d[u, u, u_I] as a yellow cross.
            plt.plot(control_point[0], control_point[1], 'X', color = 'r', markersize = 8) #Highlight the control point as a red cross.
        
        if interpolate:
            ip0,ip1=zip(*self.data_points)
            plt.plot(ip0, ip1, 'go', color = 'g')
        
        plt.show()
        
    def spline_from_basis_func(self, n):
        """
            This function evaluates S(u) for each u using the basis functions 
            given from basis_func
        """

        #Create empty list for storing L basis functions
        Ni = (len(self.u_vector)-2)*[None]
        
        #Create the basis functions for each value j (each possible hot interval) and store them in a list
        for j in range(len(self.u_vector)-2):
            Ni[j] = self.basis_func(j)
            
        #Use the control points and the basis functions to create s(u)
        
        self.u = np.linspace(min(self.u_vector)+0.001,max(self.u_vector)-0.001,n)  #Generate n-vector of u-values to generate spline
        Spline = np.empty([np.shape(self.d_vector)[1],n]) #Create empty vector for filling with points on the spline
        
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
    
    def generate_blossoms(self, u_blossom):
        i = int(self.u_vector.searchsorted(u_blossom)) #Find the "hot interval" for u_blossom
        uu = self.u_vector[i-3:i+3] #Extract the relevant interval of knot points  
        dd = self.d_vector[i-3:i+1] #Extract the relevant interval of control points
        blossoms, control_point, d43 = self.deBoor(dd,uu, u_blossom, get_blossoms = True) #Generate the blossoms and control point for u_blossom
        return blossoms, control_point, d43
    
    def load_default_data(self,mode = None):
        if mode=='d_vector':
            return np.array([[-12.73564, 9.03455],
                            [-26.77725, 15.89208],
                            [-42.12487, 20.57261],
                            [-15.34799, 4.57169],
                            [-31.72987, 6.85753],
                            [-49.14568, 6.85754],
                            [-38.09753, -1e-05],
                            [-67.92234, -11.10268],
                            [-89.47453, -33.30804],
                            [-21.44344, -22.31416],
                            [-32.16513, -53.33632],
                            [-32.16511, -93.06657],
                            [-2e-05, -39.83887],
                            [10.72167, -70.86103],
                            [32.16511, -93.06658],
                            [21.55219, -22.31397],
                            [51.377, -33.47106],
                            [89.47453, -33.47131],
                            [15.89191, 0.00025],
                            [30.9676, 1.95954],
                            [45.22709, 5.87789],
                            [14.36797, 3.91883],
                            [27.59321, 9.68786],
                            [39.67575, 17.30712]])
        if mode=='u_vector':
            KNOTS = np.linspace(0,1,26)
            KNOTS[ 1] = KNOTS[ 2] = KNOTS[ 0]
            KNOTS[-3] = KNOTS[-2] = KNOTS[-1]
            
#            print(KNOTS)            KNOTS[3:7] = 0.15

#            KNOTS[3:7] = 0.15
#            KNOTS[7:11] = 0.30
#            KNOTS[11:15] = 0.45
#            KNOTS[15:19] = 0.60            
            
#            KNOTS[3:6] = 0.15
#            KNOTS[6:9] = 0.25
#            KNOTS[9:12] = 0.35
#            KNOTS[12:15] = 0.45

            
            
            return KNOTS
        
