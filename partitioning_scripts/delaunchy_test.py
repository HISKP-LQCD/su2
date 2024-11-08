### calculates the lookup tables for the neighrest neighbor partitionings ###

import numpy as np
import HLGTTools.operators as ho
#import HLGTTools.operators.DJT.S3_sphere.partition as S3_partition
import pandas as pd
import argparse
import scipy.spatial 
### read in given arguments for script ###
parser = argparse.ArgumentParser(prog = "calculate_tables", description="""calulate the lookup tables for the partitionings""")
parser.add_argument('-m', '--m',type=int,  help = "The m/N argument for the calculation of the partitoning")
parser.add_argument('-partitioning', '--wanted_partitioning', type=str, help = "The wanted partioning")
args = parser.parse_args()
wanted_partitioning = args.wanted_partitioning
m = args.m

### calculate partitioning and weights ###
if wanted_partitioning == "Genz":    
    partitioning = ho.getSu2GenzPartitioning(m)
elif wanted_partitioning == "linear":
    partitioning = ho.getSu2LinearPartitioning(m)
elif wanted_partitioning == "Fibonacci":
    partitioning = ho.getSU2FibonacciPartitioning(N = m)
elif wanted_partitioning == "Rsc":
    partitioning = ho.getSu2RscPartitioning(N = m)
elif wanted_partitioning == "Rfcc":
    partitioning = ho.getSu2RfccPartitioning(N = m)
elif wanted_partitioning == "Volleyball":
    partitioning = ho.getSu2VolleyballPartitioning(m = m)
else:
    raise Exception("paritioning not implimented")

print("generated partitionings")

### uses Timos code to get the weights ###
weights = ho.getSU2TriangulatedIntegrationWeights(points=partitioning)

### setup the Delaunay triangulation ###
Delaunchy = scipy.spatial.Delaunay(points = partitioning)
helpresult = Delaunchy.vertex_neighbor_vertices
helpresult0 = np.copy(helpresult[0])
helpresult1 = np.copy(helpresult[1])
Delaunchy.close()
print("generated neighrest neighbors")

### create file and write the table ###
with open('lookuptable_nn.csv', 'x') as file: # note: produces an error, if the file already exists
    counter = 0 # gives the line - 1
    while counter < len(partitioning): # writes all partitionings into the file
        file.write(str(counter)) # writes the line number
        file.write(',') # seperation charater
        for element in (partitioning[counter, :]): # writes all partion coordinates
            file.write(str(element))
            file.write(",")
        file.write(str(weights[counter])) # writes the weights
        file.write(',')
        neighborarray = helpresult1[helpresult0[counter]:helpresult0[counter+1]] # creates an array with indeces of neighrest neighbors
        # loops through length of partitioning and writhe all indeces
        anothercounter = 0 
        while (anothercounter < len(weights)):
            if (anothercounter < len(neighborarray)):
                file.write(str(neighborarray[anothercounter]))
            file.write(',') # note: write this out to have a standartised file length
            anothercounter += 1
        file.write('\n') # begin next line
        counter += 1
    file.close()
