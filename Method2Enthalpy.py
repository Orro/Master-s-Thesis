# Script that calculates the potentials (Coulomb and L-J)
#
# between a water molecule and the rest of molecules of a 
#
# protein. The average potentials are calculated by dividing
#
# the totals between the number of frames.
#
########################################################

from __future__ import division
from time import time
t0 = time()

cmd.delete ("all") # Clear screen

# Variables that hold the value of the elec.and vdW potentials for every round on a loop
coulomb = 0
VLJ = 0

NumFra = 1000 # Variable that selects the number of frames

cmd.load ("cluster.pdb") # Load clusters

# Arrays to store the residues surrounded by water
stored.residues = []
stored.atoms = []
stored.numbers = []
stored.IDs = []

# Hashmaps for vdw parameters (A,B) and atom types (C)
A = {}
with open('Qoplsaa.vdw.prm', 'r') as f:
    for line in f:
        slineA = line.split()
        A[slineA[0]] = slineA[1]

B = {}
with open('Qoplsaa.vdw.prm', 'r') as f:
    for line in f:
    	slineB = line.split()
        B[slineB[0]] = slineB[2]

C = {}
with open('Qoplsaa.vdw', 'r') as f:
    for line in f:
    	slineC = line.split()
        C[(slineC[0], slineC[1])] = slineC[2]

D = {}
with open('Qoplsaa.txt', 'r') as f:
    for line in f:
    	slineD = line.split()
        D[(slineD[0], slineD[1])] = slineD[2]

# Main loop: load files, perform calculations
for i in range(1,NumFra+1):
    t2 = time()
    cmd.load("frame"+str(i)+".pdb")
    cmd.iterate("(resn HOH & name O & frame"+str(i)+") within 1 of 	cluster","stored.numbers.append(resi)") # Store all O atoms in a list
    
    # Remove membrane
    cmd.remove("resn POP")
    cmd.remove("!(byres(all within 20 of cluster))")
    
    if stored.numbers != []:
        cmd.select("totalcluster"+str(CountClusMem), "resi"+stored.numbers[0] + " & 	frame"+str(i))
        if CountClusMem > 0:
            cmd.load("totalcluster.pdb")
        # Saves the frame number
        cmd.alter("totalcluster"+str(CountClusMem), "segi = i")
	
        # Saves the cluster into a pdb file
        cmd.save("totalcluster.pdb", "totalcluster"+str(CountClusMem) +" or  		totalcluster")
        cmd.remove("totalcluster")
        cmd.delete("totalcluster")
            
        cmd.iterate("all & !(cluster) & !(resi 	"+stored.numbers[0]+")","stored.IDs.append(ID)")
	
        for j in stored.IDs:
            RO = cmd.distance("resi "+stored.numbers[0]+" & name O & frame"+str(i), "id               
             "+str(j))
            RH1 = cmd.distance("resi "+stored.numbers[0]+" & name H1 		& frame"+str(i), "id "+str(j))
            RH2 = cmd.distance("resi "+stored.numbers[0]+" & name H2 		& frame"+str(i), "id "+str(j))
            if (RO > 0) & (RH1 > 0) & (RH2 > 0):
               cmd.iterate("id "+str(j),"stored.residues=resn")
               cmd.iterate("id " +str(j),"stored.atoms=name")
	
               # vdW potential
               AOt = float(A[C[("HOH","O")]])*float(A[C[(stored.residues,stored.atoms)]])
               BOt = float(B[C[("HOH","O")]])*float(B[C[(stored.residues, stored.atoms)]])
               RO6 = RO**6
               VLJ = VLJ + AOt/RO6**2 - BOt/RO6

               AH1t = float(A[C[("HOH","H1")]])*float(A[C[(stored.residues,stored.atoms)]])
               BH1t = float(B[C[("HOH","H1")]])*float(B[C[(stored.residues, stored.atoms)]])
               RH16 = RH1**6
               VLJ = VLJ + AH1t/RH16**2 - BH1t/RH16

               AH2t = float(A[C[("HOH","H2")]])*float(A[C[(stored.residues,stored.atoms)]])
               BH2t = float(B[C[("HOH","H2")]])*float(B[C[(stored.residues, stored.atoms)]])
               RH26 = RH2**6
               VLJ = VLJ + AH2t/RH26**2 - BH2t/RH26

               # Electrostatic potential
               O = float(D[("HOH","O")]) * float(D[stored.residues,stored.atoms])/RO
               H1 = float(D[("HOH","H1")]) * float(D[stored.residues,stored.atoms])/RH1
               H2 = float(D[("HOH","H2")]) * float(D[stored.residues,stored.atoms])/RH2
               coulomb = coulomb + O
               coulomb = coulomb + H1
               coulomb = coulomb + H2
				
               cmd.delete("dist*") # Delete distancies from PyMOL memory
        
    stored.numbers = []
    stored.IDs = []
    stored.atoms = []
    stored.residues = []
			
    cmd.remove("frame"+str(i)) # Remove frames from PyMOL memory
    cmd.delete("frame"+str(i))

print "Total Coulomb's potential: ", 332*coulomb, "kcal/mol"
print "Total Lennard-Jones potential: ", VLJ, "kcal/mol"

print "Coulomb's potential for my cluster: ", 332*coulomb/NumFra, "kcal/mol"
print "Lennard-Jones potential for my cluster: ", VLJ/NumFra, "kcal/mol"

t1 = time()
print "Script running time:", t1-t0
