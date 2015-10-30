# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:00:22 2015

@author: Claire
"""
import os
from copy import copy
import pandas as pd

#%%# Wrapper functions for metabolomics pre processing steps
def pick_peaks(seq_df, mode, params, data_directory, working_directory):
    # Wrapper for call to R code that uses xcms() to pick peaks in all files
    # fnames is a list with the names of all files to be analyzed
    # params is a dictionary of parameters to give to xcsm() code
    # working_directory and mode creates the rimage_file name the rimage to be saved
    # data_directory is where the mzML files are located
    # seq_df has mode and file name for each sample

    ppm = str(params['ppm'])
    snthresh = str(params['snthresh'])
    prefilter_min = str(params['prefilter_min'])
    prefilter_max = str(params['prefilter_max'])
    integrate = str(params['integrate'])
    peakwidth_min = str(params['peakwidth_min'])
    peakwidth_max = str(params['peakwidth_max'])
    noise = str(params['noise'])
        
    # Extract all the right file names
    if mode == 'negative':
        files = seq_df[seq_df['ionmode'] == 'negative']['file name']
        files = [os.path.join(data_directory, f) + '.threshold1000.mzML' for f in files]
    elif mode == 'positive':
        files = seq_df[seq_df['ionmode'] == 'positive']['file name']
        files = [os.path.join(data_directory, f) + '.mzML' for f in files]
    else:
        raise NameError('No mode specified. Cannnot process files.')

    
    # Write the file names to a file so that R can read it in
    mzdatafiles = os.path.join(working_directory, 'mzdatafiles.' + mode + '.txt')
    with open(mzdatafiles, 'w') as f:
        f.write('fnames\n')
        f.write('\n'.join(files))
    
    # Write the sample IDs to a file so that R can read it in and label the xcmsSet object with them
    sids = seq_df[seq_df['ionmode'] == mode].index
    sidsfile = os.path.join(working_directory, 'sampleIDs.' + mode + '.txt')
    with open(sidsfile, 'w') as f:
        f.write('sids\n')
        f.write('\n'.join(sids))
    
    # Define the rimage and pdf files to save results of peak-picking
    rimage_file = os.path.join(working_directory, working_directory.split('/')[-1] + '.picked_peaks.' + mode + '.Rimage')
    pdf_file = os.path.join(working_directory,working_directory.split('/')[-1] + '.picked_peaks.' + mode + '.pdf')
    
    # Call R script that's a wrapper to xcsm() peak picking
    # This R file uses xcmsSet to pick peaks from the mzML files specified in fnames.
    # It saves the output diagnostics as a pdf file, and the resulting xs xcmsSet object in the rimage file.
    # It also labels each file in the xs object by its sampleID (using the sampclass method)    
    os.system('Rscript pick_peaks.R -p ' + ppm + ' -s ' + snthresh + ' --filterMin ' + prefilter_min +
              ' --filterMax ' + prefilter_max  + ' -i ' + integrate + ' --peakMin ' + peakwidth_min + 
              ' --peakMax ' + peakwidth_max + ' -n ' + noise + ' -f ' + mzdatafiles + 
              ' --pdf ' + pdf_file + ' --rimage ' + rimage_file + ' --sids ' + sidsfile)
    
    # Track processing that's been done on samples by writing a new sequence file for just these samples
    tmp_seq_df = copy(seq_df[seq_df['ionmode'] == mode])
    tmp_seq_df['mzML files'] = files  
    tmp_seq_df['pick_peaks_rimage'] = len(files)*[rimage_file]
    print(tmp_seq_df)
    
    proc_tracker = os.path.join(working_directory, working_directory.split('/')[-1] + '.processing_tracker.' + mode + '.txt')
    tmp_seq_df.to_csv(proc_tracker, sep='\t')
    
    return rimage_file, proc_tracker

def align_peaks(rimage, batch, samples, mode, proc_file, working_directory):
    # Wrapper to R script that aligns peaks, fills peaks, and finds adducts and isotopes
    # rimage is the full path to the Rimage file that contains the xcmsSet object with all the picked peaks
    # rimage file should be for the right mode (duh), so the xs object in it should have all of the mzML files from the respective pick_peaks call
    # batch is the name of the batch to process
    # samples is a list of samples in that batch. Should be the sample IDs.
    # proc_file is a copy of the sequence file for just the files in the respective rimage file. it has extra columns indicating whether peaks have been picked, and what batches have been finished

    ## Define output file names
    all_peaks = os.path.join(working_directory, working_directory.split('/')[-1] + '.all_peaks.batch_' + batch + '.' + mode + '.csv')
    aligned_table = os.path.join(working_directory, working_directory.split('/')[-1] + '.aligned_table.batch_' + batch + '.' + mode + '.csv')

    ## Read in the proc_file
    proc_df = pd.read_csv(proc_file, sep='\t', index_col=0)
    proc_df['batch ' + batch] = [1 if s in samples else 0 for s in proc_df.index]
    
    ## Write in a temp file indicating 'batch' column for each sample as yes or no
    # The xcmsSet object in R has each sample labeled by its sampleID. The sampleIDs are in the same order as in proc_df
    tmp_file = os.path.join(working_directory, 'tmp_batch_index.txt')
    with open(tmp_file, 'w') as f:
        f.write('sample\tinbatch\n')
        for sid in proc_df.index:
            f.write(sid + '\t' + str(proc_df.loc[sid, 'batch ' + batch]) + '\n')            
    # The R script reads in this tmp_batch_index.txt file to get the classlist
    # This file indicates (with 1's and 0's) which samples to align

    os.system('Rscript align_peaks.R --rimage ' + rimage + ' --batch ' + tmp_file + ' --mode ' + mode + 
              ' --aligned ' + aligned_table + ' --allpeaks ' + all_peaks)
    
    # Update the processing tracker file
    proc_df['batch ' + batch] = [aligned_table if s in samples else 0 for s in proc_df.index]
    proc_df.to_csv(proc_file, sep='\t')
    
    ## TODO: open the aligned feature table and replace the sample names (which are currently hte mzML files) with their sample ID (???)
    aligned_df = pd.read_csv(aligned_table, sep=',')
    new_cols = list(aligned_df.columns)
    new_cols[9:9+len(samples)] = samples
    aligned_df.columns = new_cols
    aligned_df.to_csv(aligned_table + 'test', sep=',')


#%%# Other helpful functions for metabolomics preprocessing
def extract_batches(seq_df, mode):
    # Extract the batches of samples we want to align.
    # Batches should be specified in the sequence file in the column 'batches'
    # If one sample should be in multiple batches, the multiple batches should be comma-separated
    # Returns batches, a list of the batch names to process that have the right mode
    #         b2s, a dictionary with {batch: [samples in that batch]}. Samples are labeled by their sample ID in the sequence file (which are the indices of seq_df)

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
    
    ## Create dictionary with samples in each batch b2s[batch] = [samples in that batch]
    b2s = {}
    for batch in batches:
        b2s[batch] = [s for s in s2b if batch in s2b[s]]
        
    ## Extract the batches that have all samples with the specified mode
    b2m = {}
    for batch in batches:
        b2m[batch] = [seq_df.loc[s, 'ionmode'] for s in s2b if batch in s2b[s]]
        mode_match = sum([i == mode for i in b2m[batch]])/len(b2m[batch])   # this gets the percent of samples in each batch that match the input mode. Bc these are both ints, it's either 1 or 0.
        if mode_match != 1:
            batches.remove(batch)
    
    return batches, b2s
