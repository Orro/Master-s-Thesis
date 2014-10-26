# Script that obtains rotation matrices after alignment of a reference (cluster center) 
#
# and the rest of the molecules in the cluster. Then the Euler angles are calculated,
#
# stored in a list and used to obtain a histogram whose area is calculated by the trap
#
# method.
#
#######################################################################################

import sys
from time import time
import math
from numpy import *
import numpy as numpy


# Calculate Euler angles of a cluster

# Lists that store the values of psi, phi and theta
philist1 = []
philist2 = []

thetalist1 = []
thetalist2 = []

psilist1 = []
psilist2 = []

# File with the center of cluster used as reference to obtain a rotation matrix
cmd.load("cluster.pdb")

# File with the molecules of the cluster
cmd.load("totalcluster.pdb")

f3 = open( "totalcluster.pdb", "r") # Opening the text file
line3 = f3.readline() # Reading the whole line
while line3: # Loop through the txt file, line after line
    sline3 = line3.split() # Breaks the string into a list, separating it at " "
	
    atom = sline3[0]
	
    if atom == 'ATOM':

        type = sline3[2]

        if type == 'O':

            atomid = sline3[1]
		
            # Variables that hold the id for H1 and H2	
            atomidH1 = int(atomid) + 1
            atomidH2 = int(atomid) + 2
			
            cmd.select("beta", "(id "+str(atomid)+" + id "+ str(atomidH1)+" + id "+ str(atomidH2)+") and totalcluster")
            cmd.align("cluster","beta")
            trans = cmd.get_object_matrix("cluster") # Get rotation matrices

            M = [[trans[0],trans[1],trans[2]],[trans[4],trans[5],trans[6]],[trans[8],trans[9],trans[10]]]
            
            if (M[2][0] != 1 and M[2][0] != -1):
              theta1 = -math.asin(M[2][0])
              theta2 = math.pi - theta1
              psi1 = math.atan2(M[2][1]/math.cos(theta1),M[2][2]/math.cos(theta1))
              psi2 = math.atan2(M[2][1]/math.cos(theta2),M[2][2]/math.cos(theta2))
              phi1 = math.atan2(M[1][0]/math.cos(theta1),M[0][0]/math.cos(theta1))
              phi2 = math.atan2(M[1][0]/math.cos(theta2),M[0][0]/math.cos(theta2))
	
	# Store angle values on lists					
              philist1.append(phi1)
              philist2.append(phi2)
              thetalist1.append(theta1)
              thetalist2.append(theta2)
              psilist1.append(psi1)
              psilist2.append(psi2)
					
            else:
             
              phi = 0
              
              if (M[2][0] == -1):
                theta = math.pi/2
                psi = phi + math.atan2(M[0][1],M[0][2])
                philist1.append(phi)
                philist2.append(phi)
                thetalist1.append(theta)
                thetalist2.append(theta)
                psilist1.append(psi)
                psilist2.append(psi)
					
              else:
                theta = -math.pi/2			
                psi = -phi + math.atan2(-M[0][1],-M[0][2])
                philist1.append(phi)				
                philist2.append(phi)
                thetalist1.append(theta)
                thetalist2.append(theta)
                psilist1.append(psi) 
                psilist2.append(psi)

    line3 = f3.readline()

f3.close()

def Histogram(listElements, typeOfangle, set):
    '''
    Function Histogram (calculates a normed histogram and its area with the trapezoidal rule)
    Select typeOfangle: 1 = phi, 2 = theta, 3 = psi. Set selects the group of angles: 1 or 2.
    '''
    
    arrayHist = array(listElements)
    
    hist, bin_edges = numpy.histogram(arrayHist, normed = True)

    widths = diff(bin_edges)
    
    Area = (hist * widths).sum()
    Trap = numpy.trapz(hist, dx = widths[0])

    if typeOfangle != 2:
        # Call function Pxlnxpx to calculate area under the curve f(x) = -p(x)ln(p(x))
        Pxlnxpx(hist, widths)
    else:
        PxlnxpxTheta(hist, bin_edges, set)


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

def PxlnxpxTheta(histPx, bin_edgesPx, set):
    '''
    Function PxlnxpxTheta (calculates the area under the curve f(x) = p(x)ln(p(x)sin(x))
    '''
    
    pPrime = []
    theta = []
    AreaPrimex = []

    for i in range(0, len(histPx)):
        
        a = bin_edgesPx[i]
        b = bin_edgesPx[i+1]
        
        if set == 1:
            theta.append(((a+b)/2)+(math.pi/2))
            pPrime.append(histPx[i]/math.sin(theta[i]))
        else:
            theta.append(((a+b)/2)-(math.pi/2))
            pPrime.append(histPx[i]/math.sin(theta[i]))
    
    for i in range(0, len(pPrime)):
    
        if pPrime[i] > 0:
            
            AreaPrimex.append((float(pPrime[i]) * math.log(float(pPrime[i])) *    	math.sin(math.sin(theta[i]))))
                          
    widthspP = diff(bin_edgesPx)
                
    TrapPx = numpy.trapz(AreaPrimex, dx = widthspP[0])
                
    print "Area p(x)ln(p(x))sin(x): ", TrapPx

#############################################################################################

print "\n\n"
print "#### Areas for phi1, theta1 and psi1 ####"
Histogram(philist1, 1, 1)
Histogram(thetalist1, 2, 1)
Histogram(psilist1, 3, 1)

print "\n\n"
print "#### Areas for phi2, theta2 and psi2 ####"
Histogram(philist2, 1, 2)
Histogram(thetalist2, 2, 2)
Histogram(psilist2, 3, 2)
