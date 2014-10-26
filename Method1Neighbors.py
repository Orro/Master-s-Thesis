# Loop that goes through all the frames and stores the water molecules
for i in range(1,NumFra):
    # Stores water molecules according to their residual identifier             
    cmd.iterate("resn HOH & frame"+str(i),"stored.numbers.append(resi)")
	 
    for j in stored.numbers:
        k = cmd.select("resn HOH within 1 of (resi "+j + " & frame"+str(i)+")")	
		
        if k > max:
	cluster_waters.append([j,i,k])
	# Remove elements in stored.numbers
	stored.numbers = []
		
# Method that orders cluster.waters regarding the number of neighbors of each water molecule
cluster_waters_sorted = sorted(cluster_waters, key=lambda tup: tup[2], reverse = True)

counter = 0

# Loop that goes through the list of clusters and saves those with # neighbours over cut-off
while len(cluster_waters_sorted) > 0:

	d = cluster_waters_sorted[0][0]
	e = cluster_waters_sorted[0][1]
	cmd.select("cluster"+str(counter), "resi "+d + " & frame"+str(e))
	if counter > 0:
	      cmd.load("cluster.pdb")
	
	# Saves the % of conservation in the b-factor column
	cons_ratio = (cluster_waters_sorted[0][2]/NumFra)
	cmd.alter("cluster"+str(counter), "b = cons_ratio")
	cmd.alter("cluster"+str(counter), "segi = e")
	
	# Saves the cluster into a pdb file
	cmd.save("cluster.pdb", "cluster"+str(counter) +" or 	cluster")	
	cmd.remove("cluster")
	cmd.delete("cluster")	
	counter = counter + 1

	# Removes first element of cluster_waters_sorted and all water molecules within 	1 A surrounding it
	cmd.select("selecCluster", "resn HOH within 1 of resi  	"+str(cluster_waters_sorted[0][0]))
	cmd.remove("selecCluster")

	# Loop that checks if all elements of cluster_waters have the highest number of 	neighbors
	for p in range(1,len(cluster_waters_sorted)):
                  m = cluster_waters_sorted[p][0]
	     i = cluster_waters_sorted[p][1]
	     l = cmd.select("resn HOH within 1 of (resi "+m +" & frame"+str(i)+")")	     cmd.select("cluster"+str(p), "resi "+m + " & frame"+str(i))
		
	     if l > max:
	         max_neigh.append(cluster_waters_sorted[p])
	
	cluster_waters_sorted = []      # Removes elements in cluster_waters_sorted
		
	max_neigh_sorted = sorted(max_neigh, key=lambda tup: 	tup[2],reverse = True)
	
	cluster_waters_sorted = max_neigh_sorted 		

	max_neigh = []
