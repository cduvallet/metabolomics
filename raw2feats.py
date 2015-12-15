"""
Created on Wed Oct 21 21:16:48 2015

@author: Claire Duvallet

raw2feats.py converts raw mzML files into aligned feature tables.

It interfaces with a summary file named summary_file.txt as parsed by SummaryParserMtab.
The SummaryParserMtab.py code should be either in the same directory as raw2feats.py or in your PYTHONPATH variable.
The summary_file.txt contains tab-delimited fields, as further described in the documentation.

This code calls pick_peaks.R to pick peaks in the specified ionzation mode. 
All files with the given ionization mode in the sequence file must be in the data directory.
If ionization mode is negative, mzML files should be name <fname>.threshold1000.mzML
If ionization mode is positive, mzML files should be name <fname>.mzML
Where <fname> is the file name in the sequence file.

If an Rimage file with picked peaks is given in the summary_file, peak-picking with xcmsSet is skipped.
The Rimage file should contain an xcmsSet object called xs

Peaks are aligned across all batches specified in the sequence file

All output files are labeled with the mode and dataset_ID, and saved in the output_directory
If no output directory is specified, all file are saved in <homedirectory>/proc/<dataset_ID>_proc_mtab/
"""


#%%
import preprocessing_mtab as mtab
import pandas as pd
import os
import shutil
from optparse import OptionParser
from SummaryParserMtab import *

## Have the same format as raw2otu.py in terms of interfacing with directories, etc
usage = "%prog -i INPUT_DIR -o OUTPUT_DIR_FULL_PATH"
parser = OptionParser(usage)
parser.add_option("-i", "--input_dir", type="string", dest="input_dir")
parser.add_option("-o", "--output_dir", type="string", dest="output_dir", help='Full path to output directory')
parser.add_option("-r", "--raw_data", dest="raw_data", default='False', help='if True, raw data needs to be converted to mzML using MSConvert. If False, assumes mzML files already exist')
(options, args) = parser.parse_args()


if( not options.input_dir ):
    parser.error("No data directory specified.")
options.input_dir = os.path.normpath(options.input_dir)

# Parse summary file for dataset ID
summary_file = os.path.join(options.input_dir, 'summary_file.txt')
summary_obj = SummaryParser(summary_file)
summary_obj.ReadSummaryFile()
dataset_ID = summary_obj.datasetID

# Pipe stdour and stderr to logfiles in the new directory
#sys.stdout = open('/home/ubuntu/logs/stdout_' + dataset_ID + '_proc_mtab.log','w')
#sys.stderr = open('/home/ubuntu/logs/stderr_' + dataset_ID + '_proc_mtab.log','w')
#@ THOMAS: the following commented lines give me an error. Not sure what package it needs/what it does?
#def warning(*objs):
#    print("WARNING: ", *objs, file=sys.stderr)

# If no output directory specified, default to $home/proc/
homedir = os.path.expanduser("~")
if( not options.output_dir ):
    proc_dir = os.path.join('proc', dataset_ID + '_proc_mtab')
    options.output_dir = os.path.join(homedir, proc_dir)
    print("No output directory name specified.  Writing to " + options.output_dir)
else:
    options.output_dir = os.path.normpath(options.output_dir)
    print("Saving to output directory " + options.output_dir)
    
# Make a directory for the mtab processing results
working_directory = options.output_dir
try:
    os.system('mkdir ' + working_directory)
except:
    print("Processing directory for this dataset already exists.  Overwriting its contents.")

## Parse summmary file.
# Check for presence of sequence file - raise error if no sequence file provided
try:
    sequence_file_path = os.path.join(options.input_dir, summary_obj.attribute_value_mtab['SEQUENCE_FILE'])
except:
    raise NameError('No sequence file specified.')

# Sequence file separator indicates whether sequence is csv or tab-delimited. Defaults to csv.
try:
    sequence_file_separator = summary_obj.attribute_value_mtab['SEQUENCE_FILE_SEPARATOR']
except:
    sequence_file_separator = ','

# If a separate data directory is given, otherwise defaults to input directory
try:
    data_directory = summary_obj.attribute_value_mtab['DATA_DIRECTORY']
except:
    data_directory = options.input_dir

# Get whether input is raw data or already converted mzML files
raw_data = options.raw_data

# Get mode. Summary file returns either 'negative' or 'positive', regardless of whether input was pos/neg/positive/negative/, etc
try:
    mode = summary_obj.attribute_value_mtab['MODE']
except:
    raise NameError('No mode specified.')
if mode != 'negative' and mode != 'positive':
    raise NameError('Incorrect mode specified. Valid options are "negative" and "positive"')

# Get Rimage file if it is given. If an Rimage file is given, your sequence file should only contain the samples in this Rimage file (??? TODO ???)
# The Rimage file should have an xcmsSet object called xs that contains the picked peaks from your mzML files
try:
    rimage = summary_obj.attribute_value_mtab['RIMAGE_FILE']
    rimage = os.path.normpath(rimage)
except:
    rimage = ''
#%%#0. parse sequence file
## Sequence file path is specified in summary file
# Note: first column should be the sample ID
seq_df = pd.read_csv(sequence_file_path, sep=sequence_file_separator, index_col=0)
seq_df.columns = [str(col).lower() for col in seq_df.columns]

# Check that the required columns are in sequence file.
required_cols = ['file name', 'ion mode', 'batches']
for col in required_cols:
    if col not in seq_df.columns:
        raise NameError('Required column ' + col + ' not found in sequence file')

#1. MSConvert (if raw data given)
# Not supported
#if raw_data == 'True':
#    print('Need to convert to mzML files. Running MSConvert')
#    # Run MSConvert

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
    rimage, proc_file = mtab.pick_peaks(seq_df, mode, params, data_directory, working_directory, dataset_ID)
    print('[[Peak Picking]] Picking peaks in ' + mode + ' mode complete.')
    summary_obj.attribute_value_mtab['RIMAGE_FILE'] = rimage
else:
    print('[[Peak picking]] Rimage ' + rimage + ' specified. Continuing with alignment')
    # Copy the rimage file into the output directory, it it's not already there
    os.system('copy ' + rimage + ' ' + working_directory)
    # Create a proc_file from the seq_df
    proc_file = os.path.join(working_directory, dataset_ID + '.processing_tracker.' + mode + '.txt')
    proc_df = seq_df[seq_df['ion mode'] == mode]
    proc_df['rimage'] = rimage
    proc_df.to_csv(proc_file, sep='\t')

## 4. Align peaks (per batch)
# Read in the batches to be aligned from the sequence file
# batches are comma-separated in the "batches" column of seq_df
batches, b2s = mtab.extract_batches(seq_df, mode)

## Align peaks in each batch
# This R code aligns peaks, fills peaks, and finds adducts and isotopes
# It writes a aligned_peaks file and a all_peaks file, and updates the processing_file with the processed samples
for batch in batches:
    samples = b2s[batch]
    print('[[Align peaks]] Aligning batch ' + batch + ', containing samples ' + ','.join(samples)+ '...')
    mtab.align_peaks(rimage, batch, samples, mode, proc_file, working_directory, dataset_ID)
    print('[[Align peaks]] Aligning batch ' + batch + ', containing samples ' + ','.join(samples)+ '. Complete.')
    
summary_obj.attribute_value_mtab['PROCESSED'] = 'True'
summary_obj.WriteSummaryFile()

# Copy the new summary file into the output directory
shutil.copy(summary_file, os.path.join(working_directory, 'summary_file.txt'))
print('[[Processing mzML files]] Done.')