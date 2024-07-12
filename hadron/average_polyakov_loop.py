import pandas as pd
import numpy as np
import os
import yaml
import argparse
import glob

parser = argparse.ArgumentParser(
    prog='average_polyakov_loop',
    description="""format average_polyakov_loop files in the "hadron" format""")

parser.add_argument('-f', '--inputfile',
                    help="Path to the same yaml input file used for the run")
parser.add_argument('-o', '--outputfolder', help="optional output folder to save output file to")
parser.add_argument("-append_gaugemass",type=bool, default=False, help="if true, the gaugemassis appended to the output file name")
args = parser.parse_args()

# loading and parsing the input file for the needed info
nd = yaml.load(open(str(args.inputfile)), Loader=yaml.Loader)
added_string = ""
if args.append_gaugemass:
    gaugemass = nd["metropolis"]["gaugemass"]
    added_string = f'{gaugemass:.2f}'
T = nd["geometry"]["T"]
nd_omeas = nd["omeas"]
resdir = nd_omeas["res_dir"]
nd_plaquette = nd_omeas["spatial_polyakov"]
subdir = nd_plaquette["subdir"]

if args.outputfolder == None:
    outputfolder = resdir + "" + subdir + ""
else:
    outputfolder = str(args.outputfolder) + ""

def hadronize(name):
    d1 = resdir
    d2 = d1 + "" + subdir + ""
    print(d2)
    # list all files in d2 and extract configuration numbers from names
    file_names = list(glob.glob(d2+"/"+name+"_*"))
    file_index = [int(x.rsplit(name+"_", 1)[1]) for x in file_names]
    #
    # load .dat file with indices of formatted configurations, if existing
    iconf_file = d2 + "/iconfs.dat"
    # already formatted configurations
    file_index_old = []
    if os.path.exists(iconf_file):
        print("Old configurations found")
        file_index_old = pd.read_csv(
            iconf_file, sep=" ", header=None)[0].tolist()
    ####
    n_conf_old = len(file_index_old)
    print("Reading old data if present")
    df_new = pd.DataFrame(np.zeros(shape=(1, 3)),
                          columns=["i", "polyreal", "polyimag"], index=["remove"])
    df_new = df_new.astype({"i": int})
    if n_conf_old > 0:
        path_old = d2+"/"+name+".dat"
        if os.path.exists(iconf_file) and (not os.path.exists(path_old)):
            print("Error! Something bad happened, please redo the online measurements:\n",
                  iconf_file, "\nexists, but\n", path_old, "\ndoesn't.")
            raise ValueError
        ####
        df_old = pd.read_csv(path_old, sep=" ", dtype={"i": int})
        df_new = pd.concat([df_new, df_old], axis=0)
    ####
    print("loop over new and updated files")
    for i_f in range(len(file_index)):
        path_i = file_names[i_f]
        print(i_f, path_i)
        idx_i = file_index[i_f]
        df_i = pd.read_csv(path_i, sep=" ")[["polyreal", "polyimag"]]
        arr1 = df_i.to_numpy()
        if idx_i in file_index_old:
            # removing old configuration
            df_new = df_new[df_new["i"] != idx_i]
        ####
        df1 = pd.DataFrame(arr1)
        df1.insert(0, "i", idx_i)
        df1.columns = ["i", "polyreal", "polyimag"]
        # update new array with the new dataT+1
        df_new = pd.concat([df_new, df1])
        ####
    ####
    print("Cleaning up the dataframe format before exporting")
    save_iconf = True
    df_new = df_new.drop("remove")
    save_df_new = (df_new.shape[0] > 0)
    save_iconf = save_iconf and save_df_new
    # saving
    if save_df_new:
        df_new = df_new.sort_values(by=['i'])
        df_new.to_csv(outputfolder + "/"+name+added_string + ".dat", sep=" ", index=False)
    print("Saving the new iconf file")
    if save_iconf:
        pd.DataFrame(sorted(file_index)).to_csv(
            outputfolder + added_string + "iconfs.dat", header=False, sep=" ", index=False)
    ####
    print("Removing data in the old format")
    for p in file_names:
        print("removing:", p)
        os.remove(p)
    ####
####


hadronize("polyakov")