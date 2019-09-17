# -*- coding: utf-8 -*-

import numpy as np
import unittest
import curve_handler as p
import matplotlib.pyplot as plt
from timeit import default_timer as timer


# Initiate CurveDesigner
cd = p.CurveDesigner()


# Test if the norm of the diference between methods is almost equal to zero
class TestIdentity(unittest.TestCase):
    def setUp(self):
        # Take same control vector as in handler
        self.CONTROL = np.array([[-12.73564, 9.03455],
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
        self.KNOTS = np.linspace(0,1,26)
        self.KNOTS[ 1] = self.KNOTS[ 2] = self.KNOTS[ 0]
        self.KNOTS[-3] = self.KNOTS[-2] = self.KNOTS[-1]
        cd(d_vector=self.CONTROL,u_vector=self.KNOTS)
        # Do the Blossom spline
        self.blossomspline = cd.generateSpline(5000)
        # Do the basis spline
        self.basisspline = cd.splineFromBasisFunc(5000)
        # Calculate the norm of the difference between the two methods
        self.dif = np.linalg.norm(self.blossomspline-self.basisspline)
    def test(self):
        self.assertAlmostEqual(self.dif,0) #test if the norm of the diference is almost equal to zero
        

# Test if it works for a random 2d array
class TestRandom2D(unittest.TestCase):
    def setUp(self):
        self.control = np.random.rand(24,2)
        self.knots = np.linspace(0,1,26)
        self.knots[2] = self.knots[1] = self.knots[0]
        self.knots[-3] = self.knots[-2] = self.knots[-1]
        cd(self.control,self.knots)
        self.bl_spline = cd.generateSpline(5000)
        self.ba_spline = cd.splineFromBasisFunc(5000)
        self.dif = np.linalg.norm(self.bl_spline-self.ba_spline)
    def test(self): 
        self.assertAlmostEqual(self.dif,0)
 

# test if it works for a random 3d aray
class TestRandom3D(unittest.TestCase):
    def setUp(self):
        self.control = np.random.rand(24,3)
        self.knots = np.linspace(0,1,26)
        self.knots[2] = self.knots[1] = self.knots[0]
        self.knots[-3] = self.knots[-2] = self.knots[-1]
        cd(self.control,self.knots)
        self.bl_spline = cd.generateSpline(5000)
        self.ba_spline = cd.splineFromBasisFunc(5000)
        self.dif = np.linalg.norm(self.bl_spline-self.ba_spline)    
    def test(self): 
        self.assertAlmostEqual(self.dif,0)
        
        
# test if interpol points are in splines
        
class TestInterpol(unittest.TestCase):
    def setUp(self):
        self.data_points = np.random.rand(24,2)
        self.knots = np.linspace(0,1,26)
        self.knots[2] = self.knots[1] = self.knots[0]
        self.knots[-3] = self.knots[-2] = self.knots[-1]
        cd(data_points=self.data_points,u_vector=self.knots)
        self.bl_spline = cd.generateSpline(5000,'interpolate')
    def test(self): 
        self.AssertIsIn(self.data_points,self.bl_spline)


""" TODO """
# Add more relevant tests, eg diferent dimensions, number u, number of d, hot interval
# Maybe add tearDown() methods to the test classes?
    
    
if __name__ =='__main__':
    unittest.main()
    



""" Test time for two methds with timeit"""
# Timeit gives the total process time over the number of iterations set by number
# timeit.repeat() could be used to instead save the whole list of times meassured and pick the fastest time


import timeit

n = 10 # Number of iterations

# Set up codes
mysetup = '''
import numpy as np
import matplotlib.pyplot as plt
import curve_handler as p
control = np.random.rand(24,2)
knots = np.linspace(0,1,26)
knots[2]=knots[1]=knots[0]
knots[-3]=knots[-2]=knots[-1]
cd = p.CurveDesigner(control,knots)
'''

# Code to be meassured
method1 = "bl_spline = cd.generateSpline(5000,'control points')"
# Timeit statement 
time1 = timeit.timeit(setup = mysetup, stmt = method1, number = n )
# Second code to be meassured
method2 = "ba_spline = cd.splineFromBasisFunc(5000)"
# Timeit statement
time2 = timeit.timeit(setup = mysetup, stmt = method2, number = n)

# Print the results
print("The mean time for the blossom method is ", time1/n, "seconds.")
print("The mean time for the basis method is ", time2/n, "seconds.")


