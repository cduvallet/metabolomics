# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 21:16:48 2015

@author: Claire
"""
## Note: code to install packages: http://rpy.sourceforge.net/rpy2/doc-dev/html/robjects_rpackages.html#installing-removing-r-packages

## MSConvert step:
# We're going to rename the files with a suffix indicating whether it's 1000 threshold or NoThreshold

## raw2feats.py is a wrapper to be called by Master.py
# for now, don't worry about what flags are read by Master vs. raw2feats.py
# in this code, identify where users inputs should go (i.e. readouts from summary file, sequence file, or defaults)


## preprocessing_metabolomics.py contains the actual functions that do the heavy lifting
# for example, align_features(), pick_peaks(), etc will be in this script
# The inputs to those scripts are parsed by raw2feats.py and given as inputs to these functions
# these functions might just be fancy calls to R, either via command line or rpy2 fancy acrobatics
#%%
import preprocessing_metab as mtab
import pandas as pd
import os

# Like in the 16S proc, working_directory reflects what is in the DATASET_ID in the summary_file.txt
working_directory = '/home/claire/metabolomics_test2'

## Parse summmary file. Get the following params:
sequence_file_path = '/home/claire/metabolomics_test/test_data/test_sequence_file2.csv'
sequence_file_separator = ','    # deafults to ',' if not specified
data_directory = '/home/claire/metabolomics_test/test_data'    # defaults to working_directory if not specified
raw_data = False
mode = 'negative'     # This can be either 'neg' or 'pos' 
                 # WHen reading summary file, if it's given as any permutation of neg/negative or pos/positive, return the lowercase version: either 'negative' or 'positive'

# If you provide an Rimage file, your sequence file should contain only the samples in this Rimage file
# The rimage file should have a xcmsSet object called xs
rimage = '/home/claire/metabolomics_test/metabolomics_test.picked_peaks.negative.Rimage' #'/home/claire/metabolomics_test/metabolomics_test.picked_peaks.neg.Rimage'      # If peak-picking has already been done, can specify the Rimage in the summary file (like you can specify an OTU table)


#%%#0. parse sequence file
## Sequence file path is specified in summary file
# Note: first column should be the sample ID
seq_df = pd.read_csv(sequence_file_path, sep=sequence_file_separator, index_col=0)
seq_df.columns = [str(col).lower() for col in seq_df.columns]

# Check that the required columns are in sequence file.
required_cols = ['file name', 'ionmode', 'batches']
for col in required_cols:
    if col not in seq_df.columns:
        raise NameError('Required column ' + col + ' not found in sequence file')

#1. MSConvert (if raw data given)
if raw_data:
    print('Need to convert to mzML files. Running MSConvert')
    # Run MSConvert
    # TODO: Thomas - what's the best way to do that? On another node? we can call MSConvert from command line on a windows node

## 3. Peak picking
# Params is a dictionary of paramters to give to xcms() code
# TODO: add this as a user-specifiable input, from parsing summary file. Below are the defaults
params = {}
if mode == 'negative':
    params['ppm'] = 2
    params['snthresh'] = 10
    params['prefilter_min'] = 5
    params['prefilter_max'] = 1000
    params['integrate'] = 2
    params['peakwidth_min'] = 20
    params['peakwidth_max'] = 60
    params['noise'] = 1000
elif mode == 'positive':
    params['ppm'] = 3
    params['snthresh'] = 10
    params['prefilter_min'] = 5
    params['prefilter_max'] = 1000
    params['integrate'] = 2
    params['peakwidth_min'] = 20
    params['peakwidth_max'] = 60
    params['noise'] = 1000

# If no Rimage is specified, do peak picking:
if not rimage:
    print('[[Peak Picking]] No Rimage specified. Picking peaks in ' + mode + ' mode...')
    rimage, proc_file = mtab.pick_peaks(seq_df, mode, params, data_directory, working_directory)
    print('[[Peak Picking]] Picking peaks in ' + mode + ' mode complete.')
else:
    # Create a proc_file from the seq_df
    proc_file = os.path.join(working_directory, working_directory.split('/')[-1] + '.processing_tracker.' + mode + '.txt')
    proc_df = seq_df[seq_df['ionmode'] == 'negative']
    proc_df['rimage'] = len(proc_df.index) * [rimage]
    proc_df.to_csv(proc_file, sep='\t')
   
## TODO: Update summary_file to add Rimage that results from this peak picking
# update.SummaryFile(rimage)

## 4. Align peaks (per batch)
# 4.1 read in the batches to be aligned from the sequence file
# batches are comma-separated on the "batches" column of seq_df
batches, b2s = mtab.extract_batches(seq_df, mode)

## Align peaks in each batch
# This R code aligns peaks, fills peaks, and finds adducts and isotopes
# It writes a aligned_peaks file and a all_peaks file, and updates the processing_file with the processed samples
for batch in batches[0:2]:
    samples = b2s[batch]
    print('[[Align peaks]] Aligning batch ' + batch + ', containing samples ' + ','.join(samples)+ '...')
    mtab.align_peaks(rimage, batch, samples, mode, proc_file, working_directory)
    print('[[Align peaks]] Aligning batch ' + batch + ', containing samples ' + ','.join(samples)+ '. Complete.')
