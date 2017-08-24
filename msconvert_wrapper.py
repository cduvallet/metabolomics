# -*- coding: utf-8 -*-
"""
Python script to run MSConvert on the 2015-11-11 untargeted samples from WHOI

Created on Wed Nov 18 17:43:37 2015

@author: Claire
"""

import os

msconvert = r'C://Program Files/ProteoWizard/ProteoWizard 3.0.10922/msconvert.exe'

# All files in raw_data_dir and all subdirectories get processed, btw
raw_data_dir = 'C://cygwin64/home/Alm Lab/duvallet/170426_untargeted_MIT_Kuwait/raw_data'
out_data_dir = 'C://cygwin64/home/Alm Lab/duvallet/170426_untargeted_MIT_Kuwait/mzML_files'

raw_files = []
# Get all the files in the raw_data_dir
# Note: this might break if you have raw data in multiple subdirectories but maybe not.
for pathdir, dirnames, filenames in os.walk(raw_data_dir):
    for f in filenames:
        # Keep only .RAW or .raw files (not .sld files)
        if f.endswith('.raw') or f.endswith('.RAW'):
            raw_files.append(os.path.join(pathdir, f))
        
# Prepare strings for passing to command line by replacing spaces with "\ "
raw_files = ['\ '.join(i.split(' ')) for i in raw_files]
out_data_dir = '\ '.join(out_data_dir.split(' '))
msconvert = '\ '.join(msconvert.split(' '))

for f in raw_files:
    print(f),
    #Run with thresholding
    outfile = f.split('/')[-1].split('.')[0] + '.threshold1000.mzML'
    cmdstr = msconvert + ' --mzML --filter "peakPicking true 1-"  --filter "threshold absolute 1000 most-intense" -o '\
    + out_data_dir + ' --outfile ' + outfile + ' ' + f
    os.system(cmdstr)
    
    # Run without thresholding
    cmdstr = msconvert + ' --mzML --filter "peakPicking true 1-3"  --filter "msLevel 1-3" -o ' + out_data_dir + ' ' + f
    os.system(cmdstr)
    print(' .')
