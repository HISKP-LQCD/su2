### calculates the lookup tables for the neighrest neighbor partitionings ###

import numpy as np
import HLGTTools.operators as ho
import HLGTTools.operators.DJT.S3_sphere.partition as S3_partition
import pandas as pd
import argparse
import scipy.spatial 
import itertools
### read in given arguments for script ###
parser = argparse.ArgumentParser(prog = "calculate_tables", description="""calulate the lookup tables for the partitionings""")
parser.add_argument('-m', '--m',type=int,  help = "The m/N argument for the calculation of the partitoning")
parser.add_argument('-partitioning', '--wanted_partitioning', type=str, help = "The wanted partioning")
parser.add_argument("-wl", '--weightsincluded', type= int, help = "should the weights be included or set to 1", default= 1)
args = parser.parse_args()
wanted_partitioning = args.wanted_partitioning
m = args.m
weights_included = args.weightsincluded
print(weights_included)
### calculate partitioning and weights ###
if wanted_partitioning == "Genz":    
    partitioning = ho.getSu2GenzPartitioning(m)
elif wanted_partitioning == "linear":
    #partitioning = ho.getSu2LinearPartitioning(m)
    allowedCoords = [i for i in range(-m, m + 1)]
    #print(allowedCoords)
    linearLattice = np.array([
        c for c in itertools.product(allowedCoords, repeat=4)
        if (abs(c[0]) + abs(c[1]) + abs(c[2]) + abs(c[3])) == m
    ])
    partitioning_unnormalized = linearLattice 
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

#print(partitioning)
print("generated partitionings")

### uses Timos code to get the weights ###
if np.logical_and(weights_included, wanted_partitioning != "linear"):
    print("this ran")
    weights = ho.getSU2TriangulatedIntegrationWeights(points=partitioning)
elif np.logical_and(weights_included, wanted_partitioning == "linear"):
    print("hello there")
    partitioning = partitioning_unnormalized / np.linalg.norm(partitioning_unnormalized, axis=1)[:, np.newaxis]
    weights = (np.sqrt(2)/np.linalg.norm(partitioning, axis=1)[:, np.newaxis])**3
    weights = weights.flatten()
    neighborlist = []
    jetanothercounter = 0
    print(len(partitioning_unnormalized))
    while jetanothercounter < len(partitioning_unnormalized):
        helplist = []
        print("another text", partitioning_unnormalized[jetanothercounter, ])
        jetanothercounter2 = 0
        while jetanothercounter2 < len(partitioning_unnormalized):
            if ((np.abs(partitioning_unnormalized[jetanothercounter, 0]) - np.abs(partitioning_unnormalized[jetanothercounter2, 0]))**2 + (np.abs(partitioning_unnormalized[jetanothercounter, 1]) - np.abs(partitioning_unnormalized[jetanothercounter2, 1]))**2 + (np.abs(partitioning_unnormalized[jetanothercounter, 2]) - np.abs(partitioning_unnormalized[jetanothercounter2, 2]))**2 + (np.abs(partitioning_unnormalized[jetanothercounter, 3]) - np.abs(partitioning_unnormalized[jetanothercounter2, 3]))**2 == 2):
                #print(partitioning_unnormalized[jetanothercounter2, ])
                helplist.append(jetanothercounter2)
            jetanothercounter2 += 1
        neighborlist.append(helplist)
        #print(helplist)
        jetanothercounter += 1
    print(neighborlist)
else:
    weights = np.ones(shape = len(partitioning))
### setup the Delaunay triangulation ###
Delaunchy = scipy.spatial.Delaunay(points = partitioning)
helpresult = Delaunchy.vertex_neighbor_vertices
helpresult0 = np.copy(helpresult[0])
helpresult1 = np.copy(helpresult[1])
#print(Delaunchy.vertex_neighbor_vertices)
Delaunchy.close()
#print(partitioning)
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
        if wanted_partitioning != "linear":
            neighborarray = helpresult1[helpresult0[counter]:helpresult0[counter+1]] # creates an array with indeces of neighrest neighbors
        else:
            neighborarray= np.array(neighborlist[counter])
            print(neighborarray)
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
