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

working_directory = r'C:\Users\Claire\Documents\GitHub\metabolomics'

## Parse summmary file. Get the following params:
sequence_file_path = r'C:\Users\Claire\Documents\GitHub\metabolomics\test_data\test_sequence_file.csv'
sequence_file_separator = ','    # deafults to ',' if not specified
data_directory = r'C:\Users\Claire\Documents\GitHub\metabolomics\test_data'    # defaults to working_directory if not specified
raw_data = False
mode = 'negative'     # This can be either 'neg' or 'pos' 
                 # WHen reading summary file, if it's given as any permutation of neg/negative or pos/positive, return the lowercase version: either 'negative' or 'positive'

rimage = ''      # If peak-picking has already been done, can specify the Rimage in the summary file (like you can specify an OTU table)


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

##2.Grab the mzML files we'll need downstream for peak-picking
# For negative mode, get the thresholded ones. For positive mode, just the normal mzML files
files = []

if mode == 'negative':
    files = seq_df[seq_df['ionmode'] == 'negative']['file name']
    files = [os.path.join(data_directory, f) + '.threshold1000.mzML' for f in files]
elif mode == 'positive':
    files = seq_df[seq_df['ionmode'] == 'positive']['file name']
    files = [os.path.join(data_directory, f) + '.mzML' for f in files]
else:
    raise NameError('No mode specified. Cannnot process files.')


## 3. Peak picking
# Params is a dictionary of paramters to give to xcms() code
# TODO: add this as a user-specifiable input, from parsing summary file. Below are the defaults
params = {}
if mode == 'neg':
    params['ppm'] = 2
    params['snthresh'] = 10
    params['prefilter_min'] = 5
    params['prefilter_max'] = 1000
    params['integrate'] = 2
    params['peakwidth_min'] = 20
    params['peakwidth_max'] = 60
    params['noise'] = 1000
elif mode == 'pos':
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
    print('[[Peak Picking]] No Rimage specified. Picking peaks.')
    # Pick peaks for one or both modes
    if mode == 'negative':
        print('[[Peak Picking]] Picking peaks in negative mode...')
        rimage = mtab.pick_peaks(files, 'neg', params, working_directory)
        print('[[Peak Picking]] Picking peaks in negative mode complete.')
    elif mode == 'positive':
        print('[[Peak Picking]] Picking peaks in positive mode...')
        rimage = mtab.pick_peaks(files, 'pos', params, working_directory)
        print('[[Peak Picking]] Picking peaks in positive mode complete.')
    
## TODO: Update summary_file to add Rimage that results from this peak picking
# update.SummaryFile(rimage)

## 4. Align peaks (per batch)

# 4.1 read in the batches to be aligned from the sequence file
# batches are comma-separated on the "batches" column of seq_df
batches, b2s = mtab.extract_batches(seq_df, mode)

## Align peaks in each batch
# This R code aligns peaks, fills peaks, and finds adducts and isotopes
# It writes a aligned_peaks file and a all_peaks file
# The python helper function should also write a batch_description file w/ samples in batch, batch mode, parameters used to align, and output file names
print('[[Align peaks]] Aligning peaks for ' + str(len(batches)) + ' batches...')

for batch in batches:
    samples = b2s[batch]
    print('[[Align peaks]] Aligning batch ' + batch + ', containing samples ' + ','.join(samples)+ '...')
    mtab.align_peaks(batch, samples, mode, sequence_file_path)
    print('[[Align peaks]] Aligning batch ' + batch + ', containing samples ' + ','.join(samples)+ '. Complete.')


## 4. Aligning peaks (per batch)	
#4. Align peaks: For each batch specified in summary file (somehow):
#	- Load in the Rimage
#	- select a subset of xs that contains only that batch of samples
#	- run that batch through retcor.obiwarp and group.density stuffs
#4.2 Find adducts and isotopes
#	- this requires some parameters that depend on the mode!
#4.3 save aligned file and all_peaks file
#- batch samples will be specified as an input to #4.
#	- maybe I will write a file for each batch that states: which samples are in that batch, what the mode is, and what the output files are called


