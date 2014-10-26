# Script that calculates and plots the translational entropy of a cluster
#
# based on variances of coordinates. The translational entropy is calculated
#
# by obtaining a histogram whose area is calculated by the trap
#
# method.
#
#######################################################################################

import sys
from time import time
import math
from numpy import *
import numpy as numpy

# Get the coordinates

# Lists that store the values of x,y and z
xlist = []
ylist = []
zlist = []

# File with the center of cluster used as reference to obtain a rotation matrix
cmd.load("cluster963.pdb")

# File with the molecules of the cluster
cmd.load("totalcluster963.pdb")


# Open cluster and read coordinates
f3 = open( "totalcluster963.pdb", "r") # Opening the text file
line3 = f3.readline() # Reading the whole line
while line3: # Loop through the txt file, line after line
    sline3 = line3.split() # Breaks the string into a list, separating it at " "
	
    atom = sline3[0]
    
    if atom == 'ATOM':
    
        atomid = sline3[2]
        xcoord = sline3[5]
        ycoord = sline3[6]
        zcoord = sline3[7]
	
        if atomid == 'O':
        
            xlist.append(float(xcoord))
            ylist.append(float(ycoord))
            zlist.append(float(zcoord))
    
   line3 = f3.readline()
f3.close()

def Histogram(listElements):
    '''
    Function Histogram (calculates a normed histogram and its area with the trapezoidal rule)
    '''
    arrayHist = array(listElements)
    
    hist, bin_edges = numpy.histogram(arrayHist, normed = True)
    widths = diff(bin_edges)
    
    Area = (hist * widths).sum()

    Trap = numpy.trapz(hist, dx = widths[0])

    # Call function Pxlnxpx to calculate area under the curve f(x) = p(x)ln(p(x))
    Pxlnxpx(hist, widths)

def Pxlnxpx(histPx,widthsPx):
    '''
        Function Pxlnxpx (calculates the area under the curve f(x) = p(x)ln(p(x))
    '''    
    AreaPx = []

    for i in range(0, len(histPx)):
    
        if histPx[i] != 0:
            AreaPx.append((float(histPx[i])*math.log(float(histPx[i]))))
        else:
            AreaPx.append(0)
    
    TrapPx = numpy.trapz(AreaPx, dx = widthsPx[0])

    print "Area p(x)ln(p(x)): ", TrapPx

############# I/O test ################################################################################

print "\n\n"
print "#### Area ####"
print "X coordinate:"
Histogram(xlist)
print "Y coordinate:"
Histogram(ylist)
print "Z coordinate:"
Histogram(zlist)
