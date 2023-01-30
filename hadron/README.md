# Analysis scripts for `hadron`

This folder contains scripts for the online (offline) measurements with the "R" library "hadron": https://github.com/HISKP-LQCD/hadron.

- file formatting: the data is formatted such that is more easily readable.
- analysis of the measurements

## Routines available

- `glueball.py`: converts the data for the interpolators and correlators of the glueball. 
  1. This is done in a separate step and not directly in the `C++` code on purpose. One may want to generate configurations and do offline measurements afterwards in parallel (multiple jobs, one for each bunch of measurements).
  2. The scripts takes care of deleting the files in the old format. For each choice of the parameters and interpolator type, only 4 files are saved in the corresponding folders:
    
      - `iconf.dat` : list of trajectory indices that have been put in `hadron` format
      
      - `pp.dat`, `pm.dat`, `mp.dat`, `mm.dat`: correlators
    
  3. The storage space saved in this way is mostly negligible, but avoids potential filesystem slowdowns due to the access to large number of files. Besides, the analysis scripts won't have to read the configurations one by one.
  4. The script uses the same input file used for the runs
