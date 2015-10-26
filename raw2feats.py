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
mode_to_process = 'both'        # This can be either 'neg','pos', or both 
                            # WHen reading summary file, if it's given as any permutation of neg/negative or pos/positive, return the lowercase version: either neg or pos
                            # If no default specified, return 'both'
rimage = ''      # If peak-picking has already been done, can specify the Rimage in the summary file (like you can specify an OTU table)
                 # If an rimage is specified, need to also specify a mode_to_process (can't be both!)

##0. parse sequence file
#	- for each sample ID, what is the:
#		- file name
#		- mode
#	- batches of interest... in the sequence file?
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
neg_files = []
pos_files = []
if mode_to_process == 'neg':
    neg_files = seq_df[seq_df['ionmode'] == 'negative']['file name']
elif mode_to_process == 'pos':
    pos_files = seq_df[seq_df['ionmode'] == 'positive']['file name']
else:
    neg_files = seq_df[seq_df['ionmode'] == 'negative']['file name']
    pos_files = seq_df[seq_df['ionmode'] == 'positive']['file name']

# Make sure file names point to full path of file
neg_files = [os.path.join(data_directory, f) + '.mzML' for f in neg_files]
pos_files = [os.path.join(data_directory, f) + '.mzML' for f in pos_files]

# Get the thresholded files for negative mode processing
neg_files = [f.rsplit('.',1)[0] + '.threshold1000.' + f.rsplit('.',1)[-1] for f in neg_files]

## 3. Peak picking
# Params is a dictionary of paramters to give to xcms() code
# TODO: add this as a user-specifiable input, from parsing summary file. Below are the defaults
neg_params = {}
neg_params['ppm'] = 2
neg_params['snthresh'] = 10
neg_params['prefilter_min'] = 5
neg_params['prefilter_max'] = 1000
neg_params['integrate'] = 2
neg_params['peakwidth_min'] = 20
neg_params['peakwidth_max'] = 60
neg_params['noise'] = 1000

pos_params = {}
pos_params['ppm'] = 3
pos_params['snthresh'] = 10
pos_params['prefilter_min'] = 5
pos_params['prefilter_max'] = 1000
pos_params['integrate'] = 2
pos_params['peakwidth_min'] = 20
pos_params['peakwidth_max'] = 60
pos_params['noise'] = 1000

# If no Rimage is specified (or if an Rimage is specified, but no correct mode is specified), do peak picking:
if (rimage and mode_to_process != 'neg' and mode_to_process != 'pos') or (not rimage):
    print('[[Peak Picking]] No Rimage specified. Picking peaks.')
    # Pick peaks for one or both modes
    if mode_to_process == 'neg':
        print('[[Peak Picking]] Picking peaks in negative mode...')
        mtab.pick_peaks(neg_files, 'neg', neg_params, working_directory)
        print('[[Peak Picking]] Picking peaks in negative mode complete.')
    elif mode_to_process == 'pos':
        print('[[Peak Picking]] Picking peaks in positive mode...')
        mtab.pick_peaks(pos_files, 'pos', pos_params, working_directory)
        print('[[Peak Picking]] Picking peaks in positive mode complete.')
    else:
        print('[[Peak Picking]] Picking peaks in negative mode...')
        mtab.pick_peaks(neg_files, 'neg', neg_params, working_directory)
        print('[[Peak Picking]] Picking peaks in negative mode complete.')
        print('[[Peak Picking]] Picking peaks in positive mode...')
        mtab.pick_peaks(pos_files, 'pos', pos_params, working_directory)
        print('[[Peak Picking]] Picking peaks in positive mode complete.')

## 4. Align peaks (per batch)

# 4.1 read in the batches from the sequence file
# batches are comma-separated on the "batches" column

## Create dictionary of batches for each sample. s2b[sample] = [batches that sample is in]
s2b = {key: value for key, value in zip(seq_df.index, seq_df['batches'])}
batches = []
for s in s2b:
    if isinstance(s2b[s], str):
        s2b[s] = [i.strip() for i in s2b[s].strip().split(',')]
    else:
        s2b[s] = [i.strip() for i in str(s2b[s]).strip().split(',')]
    batches = batches + s2b[s]
    
# nan and 0 are not valid batches, and will be skipped
batches = list(set(batches))
if '0' in batches:
    batches.remove('0')
if 'nan' in batches:
    batches.remove('nan')

print('[[Align peaks]] Aligning peaks for ' + str(len(batches)) + ' batches...')
## Create dictionary with samples in each batch b2s[batch] = [samples in that batch]
b2s = {}
for batch in batches:
    b2s[batch] = [s for s in s2b if batch in s2b[s]]
    
## Track the mode of each batch
b2m = {}
for batch in batches:
    b2m[batch] = [seq_df.loc[s, 'ionmode'] for s in s2b if batch in s2b[s]]
    if b2m[batch].count(b2m[batch][0]) != len(b2m[batch]):
        print('[[Align peaks]] Batch ' + batch + ' contains mixed negative and positive ion modes. Skipping batch.')
        del b2m[batch]
        del b2s[batch]

## Align peaks in each batch
# This R code aligns peaks, fills peaks, and finds adducts and isotopes
# It writes a aligned_peaks file and a all_peaks file
# The python helper function should also write a batch_description file w/ samples in batch, batch mode, parameters used to align, and output file names

for batch in batches:
    mode = b2m[batch]
    samples = b2s[batch]
    mtab.align_peaks(batch, samples, mode, sequence_file_path)


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


