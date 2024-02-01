
import pandas as pd
import numpy as np
import os
import yaml
import argparse
import glob
import subprocess

parser = argparse.ArgumentParser(
	prog='gradient_flow',
	description="""format gradient flow measurements "hadron" format""")

parser.add_argument(
	'-f', '--inputfile',
	help="Path to the same yaml input file used for the run"
)
# clean flowed configurations
parser.add_argument(
	'-cfc',
	action=argparse.BooleanOptionalAction,
	help="Remove the flowed configurations"
)
# clean gradient flow
parser.add_argument(
	'-cgf',
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
tmin, tmax = nd_gflow["tstart"], nd_gflow["tmax"]
eps = nd_gflow["epsilon"]
Nt = int((tmax - tmin)/(2*eps)+1)
t_flow = [tmin + eps*i for i in range(1, Nt+1)] # flow times

def hadronize():
	print("## Formatting the flowed plaquettes P, P_ss")
	d2 = resdir+"/"+subdir+"/"
	print("## Reading data from:")
	print(d2)
	print("## listing all files in there and extract configuration numbers from names")

	file_names = subprocess.run(
		"find {d2} -type f -name \"gradient_flow*\" -and -not -name \"*.conf\" ".format(d2=d2),
		shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
	).stdout.split("\n")[:-1]
	file_names_confs = subprocess.run(
		"find {d2} -type f -name \"gradient_flow*.conf\" ".format(d2=d2),
		shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
	).stdout.split("\n")[:-1]
	file_index = [int(x.split("gradient_flow.")[-1]) for x in file_names]

	if len(file_names) == 0:
		print("## Warning: No new online measurements found in ", d2)
	####

	# sorting
	file_names = [x for _,x in sorted(zip(file_index, file_names))]
	file_index = sorted(file_index)

	if len(file_names) > n_meas:
		print("### considering only the first n_meas files")
		file_names = file_names[0:n_meas]
		file_index = file_index[0:n_meas]
	####

	Ng = len(file_index)

	## t_file = d2 + "/t.dat"
	## iconf_file = d2 + "/iconfs.dat"
	P_file =  d2 + "/P.dat"
	P_ss_file =  d2 + "/P_ss.dat"

	# df_t = pd.DataFrame({'t': t_flow})
	# df_t.to_csv(t_file, index=False)

	# df_conf = pd.DataFrame({'i': file_index})
	# df_conf.to_csv(iconf_file, index=False)

	P = np.zeros(shape=(Ng, Nt))
	P_ss = np.zeros(shape=(Ng, Nt))

	for i in range(Ng):
		f = file_names[i]
		df_i = pd.read_csv(f, sep=" ", skip_blank_lines=True)
		P_i = df_i["P"].to_numpy()[0:Nt]
		P_ss_i = df_i["P_ss"].to_numpy()[0:Nt]
		P[i,:] = P_i
		P_ss[i,:] = P_ss_i
	####

	df_P = pd.DataFrame(P, index = file_index, columns=t_flow)
	df_P.to_csv(P_file, sep=" ")

	df_P_ss = pd.DataFrame(P_ss, index = file_index, columns=t_flow)
	df_P_ss.to_csv(P_ss_file, sep=" ")

	print("## Removing data in the old format")
	if args.cfc:
		for f in file_names_confs:
			print("## removing", f)
			os.remove(f)
	####
	if args.cgf:
		for p in file_names:
			print("#### removing:", p)
			os.remove(p)
		####
	####
####


print("---")
hadronize()

