# -*- coding: utf-8 -*-
"""
Python script to run MSConvert on the 2015-11-11 untargeted samples from WHOI

Created on Wed Nov 18 17:43:37 2015

@author: Claire
"""

import os

msconvert = r'C://MSConvert/msconvert.exe'

# All files in raw_data_dir and all subdirectories get processed, btw
raw_data_dir = 'C://Users/Claire/Documents/GitHub/metabolomics/UW_2015-11-11_raw_data_ClaireOnly'
out_data_dir = 'C://Users/Claire/Documents/GitHub/metabolomics/mzML_files'

raw_files = []
# Get all the files in the raw_data_dir
# Note: this might break if you have raw data in multiple subdirectories but maybe not.
for pathdir, dirnames, filenames in os.walk(raw_data_dir):
    for f in filenames:
        raw_files.append(os.path.join(pathdir, f))

for f in raw_files:
    #Run with thresholding
    cmdstr = msconvert + ' --mzML --filter "peakPicking true 1-"  --filter "threshold absolute 1000 most-intense" -o '\
    + out_data_dir + ' --outfile ' + f.split('\\')[-1].split('.')[0] + '.threshold1000.mzML ' + f
    os.system(cmdstr)
    
    # Run without thresholding
    cmdstr = msconvert + ' --mzML --filter "peakPicking true 1-3"  --filter "msLevel 1-3" -o ' + out_data_dir + ' ' + f
    os.system(cmdstr)
