# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:00:22 2015

@author: Claire
"""
import os

#%%# Wrapper functions for metabolomics pre processing steps
def pick_peaks(fnames, mode, params, working_directory):
    # Wrapper for call to R code that uses xcms() to pick peaks in all files
    # fnames is a list with the names of all files to be analyzed
    # params is a dictionary of parameters to give to xcsm() code
    # working_directory and mode creates the rimage_file name the rimage to be saved

    ppm = str(params['ppm'])
    snthresh = str(params['snthresh'])
    prefilter_min = str(params['prefilter_min'])
    prefilter_max = str(params['prefilter_max'])
    integrate = str(params['integrate'])
    peakwidth_min = str(params['peakwidth_min'])
    peakwidth_max = str(params['peakwidth_max'])
    noise = str(params['noise'])
        
    # Write the file names to a file so that R can read it in
    mzdatafiles = working_directory + '.mzdatafiles.' + mode + '.txt'
    with open(mzdatafiles, 'w') as f:
        f.write('fnames')
        f.write('\n'.join(fnames))
    
    # Define the rimage and pdf files to save results of peak-picking
    rimage_file = working_directory + '.picked_peaks.' + mode + '.Rimage'
    pdf_file = working_directory + '.picked_peaks.' + mode + '.pdf'
    
    # Call R script that's a wrapper to xcsm() peak picking
    # This R file uses xcmsSet to pick peaks from the mzML files specified in fnames.
    # It saves the output diagnostics as a pdf file, and the resulting xs xcmsSet object in the rimage file.
    # It also labels each file in the xs object by its full file name (using the sampclass method)    
    os.system('Rscript pick_peaks.R -p ' + ppm + ' -s ' + snthresh + ' --filterMin ' + prefilter_min +
              ' --filterMax ' + prefilter_max  + ' -i ' + integrate + ' --peakMin ' + peakwidth_min + 
              ' --peakMax ' + peakwidth_max + ' -n ' + noise + ' -f ' + fnames + ' --pdf ' + pdf_file + ' --rimage ' + rimage_file)
    
    return rimage_file

def align_peaks(batch, samples, mode, seq_df):
    # Wrapper to R script that aligns peaks, fills peaks, and finds adducts and isotopes
    # batch is the name of the batch to process
    # samples is a list of samples in that batch
    # seq_df is the pandas dataframe containing the sequence file. This is used to make a "class list" which xcms uses to
    # R code will need to turn the list of samples given into a classes list where all of the samples are in order **will need to check this part**

    class_list = ['smthg']
    
    os.system('Rscript path-to-r-file --flags indicating mode, class list (as a file?), ...')
    
    with open(batch + '_align_proc.txt', 'w') as f:
        f.write('Processing batch ' + batch + '\n')
        f.write('Samples:\n' + '\n'.join(samples) + '\n')
        f.write('Mode:\t' + mode + '\n')

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
