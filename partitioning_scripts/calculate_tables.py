import numpy as np
import HLGTTools.operators as ho
import HLGTTools.operators.DJT.S3_sphere.partition as S3_partition
import pandas as pd
import argparse
import scipy.spatial 

parser = argparse.ArgumentParser(prog = "calculate_tables", description="""calulate the lookup tables for the partitionings""")
parser.add_argument('-m', '--m',type=int,  help = "The m/N argument for the calculation of the partitoning")
parser.add_argument('-partitioning', '--wanted_partitioning', type=str, help = "The wanted partioning")
args = parser.parse_args()

wanted_partitioning = args.wanted_partitioning
m = args.m
print("hello there")
#wanted_partitioning = "Genz"
#m = 3
def distance_to_identity(point):
    return S3_partition.get_distance_S3(y1 = point, y2 = np.array([1, 0, 0, 0]))

def distance_to_idenity_arrays(point0, point1, point2, point3):
    helplist = []
    i = 0
    while (i < len(point0)):
        helplist.append(distance_to_identity(point=np.array([point0[i], point1[i], point2[i], point3[i]])))
        i = i +1
    return np.array(helplist)

def dagger(tree, point0, point1, point2, point3):
    helplist = []
    i = 0
    while (i < len(point0)):
        helpvalue = tree.query(np.array([point0[i], -point1[i], -point2[i], -point3[i]]))
#        print(helpvalue)
        helplist.append(helpvalue[1])
        i = i + 1
    return np.array(helplist)

def multiplication(tree, pointx, pointy):
    
    result0 = pointx[0]*pointy[0] - np.dot(pointx[1:], pointy[1:])
#    print(result0)
    result_vecotr = pointx[0] * pointy[1:] + pointy[0] * pointx[1:] - np.cross(pointx[1:], pointy[1:])
#    print(result_vecotr)
    return tree.query(np.array([result0, result_vecotr[0], result_vecotr[1], result_vecotr[2]]))[1]

def addition(tree, pointx, pointy):
    result_point = pointx + pointy
    determinant = (result_point[0]**2) + np.dot(result_point[1:], result_point[1:])
    print(result_point)
    if determinant != 0:
        result_point = result_point/np.sqrt(determinant)
    print(result_point)
    return tree.query(result_point)[1]

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



weights = ho.getSU2TriangulatedIntegrationWeights(points=partitioning)

k_d_tree = scipy.spatial.KDTree(data = partitioning)
d, i = k_d_tree.query(np.array([[0.5,0.5,0.5,-1]]))

print(multiplication(tree = k_d_tree, pointx=np.array([1,2,3,4]), pointy=np.array([1,2,3,4])))
dataframe1 = pd.DataFrame(data = partitioning)
dataframe1["weights"] = weights
dataframe1["distance_to_identitiy"] = distance_to_idenity_arrays(dataframe1[0].to_numpy(), dataframe1[1].to_numpy(), dataframe1[2].to_numpy(), dataframe1[3].to_numpy())
print(dataframe1)
dataframe1["dagger"] = dagger(tree = k_d_tree, point0 = dataframe1[0].to_numpy(), point1 = dataframe1[1].to_numpy(), point2 = dataframe1[2].to_numpy(), point3 = dataframe1[3].to_numpy())
#print(dataframe1)
dataframe1.to_csv("lookuptable1.csv", index_label="i")
#print(dataframe1)



dataframe_multiplication = pd.DataFrame()
dataframe_additon = pd.DataFrame()
indexx = 0
cooridnate0 = dataframe1[0].to_numpy()
coordinate1 = dataframe1[1].to_numpy()
coordinate2 = dataframe1[2].to_numpy()
coordinate3 = dataframe1[3].to_numpy()

while indexx < len(dataframe1):
    indexy = 0
    helplist = []
    helplist2 = []
    while indexy < len(dataframe1):
        helplist2.append(addition(tree = k_d_tree, pointx=np.array([cooridnate0[indexx], coordinate1[indexx], coordinate2[indexx], coordinate3[indexx]]), pointy = np.array([cooridnate0[indexy], coordinate1[indexy], coordinate2[indexy], coordinate3[indexy]])))
        helplist.append(multiplication(tree = k_d_tree, pointx=np.array([cooridnate0[indexx], coordinate1[indexx], coordinate2[indexx], coordinate3[indexx]]), pointy = np.array([cooridnate0[indexy], coordinate1[indexy], coordinate2[indexy], coordinate3[indexy]])))
        indexy = indexy + 1
    dataframe_multiplication[indexx] = helplist
    dataframe_additon[indexx] = helplist2
    indexx = indexx + 1
#print(addition(tree = k_d_tree, pointx = np.array([-1, 0, 0, 0]), pointy = np.array([-1, 0, 0, 0])))
print(dataframe_multiplication)
print(dataframe_additon)
dataframe_additon.to_csv("lookup_table_addition.csv", index_label = "i")
dataframe_multiplication.to_csv("lookup_table_multiplication.csv", index_label="i")