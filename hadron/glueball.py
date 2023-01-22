
import pandas as pd
import numpy as np
import os
import yaml
import argparse
import glob


parser = argparse.ArgumentParser(
    prog='glueball',
    description="""format glueball correlators and interpolators""")

parser.add_argument('-f', '--inputfile')
args = parser.parse_args()

# loading and parsing the input file for the needed info
nd = yaml.load(open(args.inputfile), Loader=yaml.Loader)

T = nd["geometry"]["T"]
nd_omeas = nd["omeas"]
resdir = nd_omeas["res_dir"]
nd_glueball = nd_omeas["glueball"]

save_corr = nd_glueball["correlator"]
nd_interpolator = nd_glueball["interpolator"]

do_APE_smearing = bool(nd_glueball["do_APE_smearing"])
n_APE = []
alpha_APE = 0.0
if (do_APE_smearing):
    nd_APE_smearing = nd_glueball["APE_smearing"]
    n_APE = nd_APE_smearing["n"]
    alpha = nd_APE_smearing["alpha"]
####

rmin = int(nd_interpolator["rmin"])
rmax = int(nd_interpolator["rmax"])
interp_type = nd_interpolator["type"]
interp_save = nd_interpolator["save"]


for i in range(rmin, rmax+1):
    d1 = resdir+"/"+"glueball/"
    d2 = d1 + "interpolator/"+interp_type+"/"
    d2 += "smearAPEn{n_APE}alpha{alpha}/".format(n_APE=n_APE, alpha=alpha)+"/"
    d2 += str(i)+"/"
    print(d2)
    # list all files in d2 and extract configuration numbers from names
    file_names = glob.glob(d2+"/phi_*")
    file_index = sorted([int(x.split("phi_")[1]) for x in file_names])
    #
    # load .dat file with indices of formatted configurations, if existing
    iconf_file = d2 + "/iconfs"
    # already formatted configurations
    file_index_old = []
    if os.path.exists(iconf_file):
        file_index_old = pd.read_csv(
            iconf_file, sep=" ", header=None)[0].tolist()
    ####
    n_conf_old = len(file_index_old)
    # new matrices for a each J^{PC}
    # dict({"pp": pd.DataFrame(), "pm": pd.DataFrame(), "mp": pd.DataFrame(), "mm":pd.DataFrame()})
    df_new = dict()
    for j_pc in ["pp", "pm", "mp", "mm"]:
        # loop over all J^{PC} quantum numbers
        df_new[j_pc] = pd.DataFrame(np.zeros(shape=(1, T+1)), columns=["i"]+[
            t for t in range(T)], index=["remove"])
        if n_conf_old > 0:
            pd.concat([df_new[j_pc], pd.read_csv(d2+"/"+j_pc, sep=" ")])
        ####
        df_new[j_pc] = df_new[j_pc].astype({"i": int})
    ####
    # loop over new and updated files
    for i_f in range(len(file_index)):
        path_i = file_names[i_f]
        idx_i = file_index[i_f]
        df_i = pd.read_csv(path_i, sep=" ")
        for j_pc in ["pp", "pm", "mp", "mm"]:
            # loop over all J^{PC} quantum numbers
            arr1 = df_i[j_pc].to_numpy()
            if idx_i in file_index_old:
                # removing old configuration
                df_new[j_pc] = df_new[j_pc][df_new[j_pc]["i"] != idx_i]
            ####
            df1 = pd.DataFrame(arr1).transpose()
            df1.insert(0, "i", idx_i)
            # update new array with the new data
            df_new[j_pc] = pd.concat([df_new[j_pc], df1])
            ####
        ####
    ####
    for j_pc in ["pp", "pm", "mp", "mm"]:
        df_new[j_pc] = df_new[j_pc].drop("remove")
        df_new[j_pc].to_csv(d2 + "/"+j_pc, sep=" ", index=False)
        for i_f in range(len(file_index)):
            # remove the old data!!!
            pass
        ####
    ####
    # save the new iconf file
    pd.DataFrame(file_index).to_csv(d2+"iconfs", header=False, sep=" ", index=False)
####


