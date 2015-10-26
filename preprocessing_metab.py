# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:00:22 2015

@author: Claire
"""
import os

## Wrapper functions for metabolomics pre processing steps

def pick_peaks(fnames, mode, params, working_directory):
    # fnames is a file with the names of all files to be analyzed
    # Wrapper for call to R code that uses xcms() to pick peaks in all files
    # params is a dictionary of paramters to give to xcsm() code
    # working directory creates the rimage_file name the rimage to be saved

    ppm = str(params['ppm'])
    snthresh = str(params['snthresh'])
    prefilter_min = str(params['prefilter_min'])
    prefilter_max = str(params['prefilter_max'])
    integrate = str(params['integrate'])
    peakwidth_min = str(params['peakwidth_min'])
    peakwidth_max = str(params['peakwidth_max'])
    noise = str(params['noise'])
        
    # Write the file names to a file so that R can read it in...
    mzdatafiles = 'mzdatafiles.' + mode + '.txt'
    with open(mzdatafiles, 'w') as f:
        f.write('fnames')
        f.write('\n'.join(fnames))
    
    # Call R script that's a wrapper to xcsm() peak picking
    os.system('Rscript pick_peaks.R -p ' + ppm + ' -s ' + snthresh + ' --filterMin ' + prefilter_min +
              ' --filterMax ' + prefilter_max  + ' -i ' + integrate + ' --peakMin ' + peakwidth_min + 
              ' --peakMax ' + peakwidth_max + ' -n ' + noise + ' -f ' + fnames + ' -w ' + working_directory)
    

def align_peaks(batch, samples, mode, sequence_file):
    # Wrapper to R script that aligns peaks, fills peaks, and finds adducts and isotopes
    # R code will need to turn the list of samples given into a classes list where all of the samples are in order **will need to check this part**

    class_list = ['smthg']
    
    os.system('Rscript path-to-r-file --flags indicating mode, class list (as a file?), ...')
    
    with open(batch + '_align_proc.txt', 'w') as f:
        f.write('Processing batch ' + batch + '\n')
        f.write('Samples:\n' + '\n'.join(samples) + '\n')
        f.write('Mode:\t' + mode + '\n')
    
