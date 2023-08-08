
import pandas as pd
import numpy as np
import os
import yaml
import argparse
import glob


parser = argparse.ArgumentParser(
    prog='gradient_flow',
    description="""format gradient flow measurements "hadron" format""")

parser.add_argument(
    '-f', '--inputfile',
    help="Path to the same yaml input file used for the run"
)
parser.add_argument(
    '-cfc', '--clean_flowed_configs',
    action=argparse.BooleanOptionalAction,
    help="Remove the flowed configurations"
)
parser.add_argument(
    '-cgf', '--clean_gradient_flow',
    action=argparse.BooleanOptionalAction,
    help="""Remove the gradient flow online measurements: 
    ACHTUNG: do it only if you don't want to flow them for larger times"""
)
parser.set_defaults(cfc=False, cgf=False)

args = parser.parse_args()

# loading and parsing the input file for the needed info
nd = yaml.load(open(str(args.inputfile)), Loader=yaml.Loader)


nd_omeas = nd["omeas"]

resdir = nd_omeas["res_dir"]

n_meas = nd_omeas["n_meas"]

nd_gflow = nd_omeas["gradient_flow"]
subdir = nd_gflow["subdir"]
tmin, tmax = nd_gflow["tmin"], nd_gflow["tmax"]
eps = nd_gflow["epsilon"]
Nt = int((tmax - tmin)/(2*eps))

def hadronize(name_obj):
    print("## Formatting the observable:", name_obj)
    d2 = resdir+"/"+subdir+"/"
    print("## Reading data from:")
    print(d2)
    print("## listing all files in there and extract configuration numbers from names")
    f1 = glob.glob(d2+"/gradient_flow"+"*")
    f2 = glob.glob(d2+"/gradient_flow*.conf")
    file_names = [f for f in f1 if not f in f2]
    file_index = [int(x.split("gradient_flow.")[1]) for x in file_names]

    if len(file_names) == 0:
        print("## Skipping: No new online measurements found in ", d2)
        return
    ####

    # sorting
    file_names = [x for _,x in sorted(zip(file_index, file_names))]
    file_index = sorted(file_index)

    if len(file_names) > n_meas:
        print("### considering only the first n_meas files")
        file_names = file_names[0:n_meas]
        file_index = file_index[0:n_meas]
    ####

    #
    iconf_file = d2 + "/iconfs.dat"
    # already formatted configurations
    file_index_old = []
    if os.path.exists(iconf_file):
        print("## load iconfs.dat file with indices of previously formatted configurations")
        print(iconf_file)
        file_index_old = pd.read_csv(
            iconf_file, sep=" ", header=None)[0].tolist()
    else:
        print("## no old configurations found")
    ##
    n_conf_old = len(file_index_old)
    # new matrices for a each J^{PC}
    df_col_names = ["i"]+[str(eps*t) for t in range(1, Nt+1)]
    df_new = pd.DataFrame(np.zeros(shape=(1, Nt+1)),
        columns=df_col_names, index=["remove"])
    df_new = df_new.astype({"i": int})
    if n_conf_old > 0:
        path_old = d2+"/"+name_obj+".dat"
        if os.path.exists(iconf_file) and (not os.path.exists(path_old)):
            print("Error! Something bad happened, please redo the online measurements:\n",
                  iconf_file, "\nexists, but\n", path_old, "\ndoesn't.")
            raise ValueError
        ####
        print("#### reading the old dataframe of the correlator")
        df_old = pd.read_csv(path_old, sep=" ", dtype={"i": int})
        df_new = pd.concat([df_new, df_old], axis=0)
    ####
    print("### loop over new/updated files")
    for i_f in range(len(file_index)):
        path_i = file_names[i_f]
        print(i_f, path_i)
        idx_i = file_index[i_f]
        df_i = pd.read_csv(path_i, sep=" ")
        arr1 = df_i.to_numpy()
        if idx_i in file_index_old:
            # removing old configuration
            df_new = df_new[df_new["i"] != idx_i]
        ####
        df1 = pd.DataFrame(arr1).transpose()
        df1.insert(0, "i", idx_i)
        df1.columns = df_col_names
        # update new array with the new data
        df_new = pd.concat([df_new, df1])
        ####
    ####
    print("### Cleaning up the dataframe format before exporting")
    save_iconf = True 
    df_new = df_new.drop("remove")
    save_df_new = (df_new.shape[0]>0)
    save_iconf = save_iconf and save_df_new
    # saving
    if save_df_new:
        df_new = df_new.sort_values(by=['i'])
        df_new.to_csv(d2 + "/"+name_obj+".dat", sep=" ", index=False)
    ####
    print("## Saving the new iconf.dat file")
    if save_iconf:
        pd.DataFrame(
            sorted(file_index_old + file_index)).to_csv(
                d2+"iconfs.dat", header=False, sep=" ", index=False)
    ####
    print("## Removing data in the old format")
    if args.cfc:
        for f in f2:
            print("## removing", f)
            # os.remove(f)
    ####
    if args.cgf:
        for p in file_names:
            print("#### removing:", p)
            #os.remove(p)
        ####
    ####
####



print("---")
hadronize("P")

print("---")
hadronize("P_ss")

