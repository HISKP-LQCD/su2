### calculates the lookup tables for the neighrest neighbor partitionings ###

import numpy as np
import HLGTTools.operators as ho
import HLGTTools.operators.DJT.S3_sphere.partition as S3_partition
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

### gives the neighrest neighbors conected via a Delaunay triagulation ###
def find_neighbors(pindex, triang):
    return triang.vertex_neighbor_vertices[1][triang.vertex_neighbor_vertices[0][pindex]:triang.vertex_neighbor_vertices[0][pindex+1]]

### uses Timos code to get the weights ###
weights = ho.getSU2TriangulatedIntegrationWeights(points=partitioning)

### setup the Delaunay triangulation ###
Delaunchy = scipy.spatial.Delaunay(points = partitioning)

### setup the dataframe that is later used to create the csv ###
dataframe1 = pd.DataFrame(data = partitioning)
dataframe1["weights"] = weights
nnarraay = np.empty(shape = (len(partitioning), len(partitioning)))
nnarraay[: ] = np.nan

### append the neighrest neighbors to the dataframe ###
counter = 0
while counter < len(partitioning):
    neighborarray = find_neighbors(counter, Delaunchy)
    nnarraay[counter, 0: len(neighborarray)] = neighborarray
    counter += 1
dataframe1 = pd.concat([dataframe1, pd.DataFrame(nnarraay)], axis = 1)

### Leaving the print command in to show the dataframe to the user before saving. This should make debugging easier ###
print(dataframe1)

### optput the dataframe as a csv ###
dataframe1.to_csv("lookuptable_nn.csv", header = False, na_rep = "")

