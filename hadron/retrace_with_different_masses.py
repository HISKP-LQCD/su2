import os
import argparse
import yaml

### load all arguments from command line ###
parser = argparse.ArgumentParser(prog = "retrace_multiple_masses", description=""" loops over retrace for multiple gauge masses""")
parser.add_argument('-f', '--inputfile', help="Path to the yaml input file for the calculations")
parser.add_argument('-numit', help="the number of loop iterations", type=int)
parser.add_argument('-startmass',default=0.0, help="the mass at which the loop starts", type=float)
parser.add_argument('-stepsize', help="the step size in the gauge masses", type=float)
parser.add_argument('-o', "--outputfolder", help = "the folder to save the outputfiles in")
parser.add_argument('-op', "--outputfolder_plaquettes", help = "The output folder for the calculated plaquettes")
parser.add_argument('-osp', "--outputfolder_polyakov", help = "the outputfolder for the spatially averaged polyakov loops")
args = parser.parse_args()

### iterate over all wanted gaugemasses ###
i = 0
while i < args.numit:
        
    #load the yaml file #
    nd = yaml.load(open(str(args.inputfile)), Loader=yaml.Loader)

    # change gauge masses #
    nd["metropolis"]["gaugemass"] = args.startmass + i*args.stepsize
    
    #save changed yaml file#
    with open(str(args.inputfile), 'w',) as f :
        yaml.dump(nd,f,sort_keys=False)
    
    #ran simulation and retrace #
    os.system("./su2-main -f " + str(args.inputfile))
    os.system(f"python3 ~/code/su2/hadron/retrace.py -f {str(args.inputfile)} --outputfolder {str(args.outputfolder)}")
    os.system(f"python3 ~/code/su2/hadron/plaquette.py -f {str(args.inputfile)} -append_gaugemass True -o {str(args.outputfolder_plaquettes)}")
    #os.system("python3 ~/code/su2/hadron/retrace.py -f " + str(args.inputfile) + "-o " + str(args.outputfolder))
    #os.system("python3 ~/code/su2/hadron/plaquette.py -f " + str(args.inputfile) + "--outputfolder" + str(args.outputfolder))
    os.system(f"python3 ~/code/su2/hadron/average_polyakov_loop.py -f {str(args.inputfile)} -o {str(args.outputfolder_polyakov)} -append_gaugemass True")
    os.system(f"mv omeas/result3p1d.u1potential.rotated.Nt16.Ns16.b3.000000.xi1.000000.nape0.alpha1.000000coarsedistance  potential_files/result3p1d.u1potential.rotated.Nt16.Ns16.b3.000000.xi1.000000.nape0.alpha1.000000coarsedistance_{str(args.startmass + i*args.stepsize)}")
    i += 1