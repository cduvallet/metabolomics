# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:00:22 2015

@author: Claire
"""
import os
from copy import copy
import pandas as pd

#%%# Wrapper functions for metabolomics pre processing steps
def pick_peaks(seq_df, mode, params, data_directory, working_directory, dataset_ID):
    """
    Picks peaks for all files in seq_df with specified mode. 
    
    Inputs:
    =======
    - seq_df:              (pandas dataframe) sequence file read into dataframe.
    - mode:                (str) ionization mode, either 'positive' or 'negative'
    - params:              (dict) dictionary containing the parameters for peak picking
    - data_directory:      (str) path to directory containing files in seq_df
    - working_directory:   (str) path to directory where rimage file will be save
    
    Returns:
    ========
    - rimage_file:    (str) full path of saved Rimage file which contains xs object with picked peaks
    - proc_tracker:   (str) full path of file containing edited seq_df with paths to files containing completed analyses

    Creates:
    ========
    - <dataset_ID>.picked_peaks.<mode>.pdf         (file) file containing pdf output from xcmsSet
    - <dataset_ID>.picked_peaks.<mode>.Rimage      (file) Rimage file contaiing xs xcmsSet object with picked peaks
    - <dataset_ID>.processing_tracker.<mode>.txt   (file) tab-delimited file with subset of seq_df with specified mode, tracks processing done

    """

    ## Get the parameters for later call to Rscript
    ppm = str(params['ppm'])
    snthresh = str(params['snthresh'])
    prefilter_min = str(params['prefilter_min'])
    prefilter_max = str(params['prefilter_max'])
    integrate = str(params['integrate'])
    peakwidth_min = str(params['peakwidth_min'])
    peakwidth_max = str(params['peakwidth_max'])
    noise = str(params['noise'])
        
    ## Extract all the right mzML file names
    # All files with the given mode that are in the sequence file are assumed to be present in the data_directory
    # Negative mode files have suffix *.threshold1000.mzML
    # Positive mode files have suffix *.mzML
    if mode == 'negative':
        files = seq_df[seq_df['ion mode'] == 'negative']['file name']
        files = [os.path.join(data_directory, f) + '.threshold1000.mzML' for f in files]
    elif mode == 'positive':
        files = seq_df[seq_df['ion mode'] == 'positive']['file name']
        files = [os.path.join(data_directory, f) + '.mzML' for f in files]
    else:
        raise NameError('No mode specified. Cannnot process files.')
    
    ## Write the file names to a file so that R can read it in
    mzdatafiles = os.path.join(working_directory, 'mzdatafiles.' + mode + '.txt')
    with open(mzdatafiles, 'w') as f:
        f.write('fnames\n')
        f.write('\n'.join(files))
    
    ## Write the sample IDs to a file so that R can read it in and label the xcmsSet object with them
    sids = seq_df[seq_df['ion mode'] == mode].index
    sidsfile = os.path.join(working_directory, 'sampleIDs.' + mode + '.txt')
    with open(sidsfile, 'w') as f:
        f.write('sids\n')
        f.write('\n'.join(sids))
    
    ## Define the rimage and pdf files to save results of peak-picking
    rimage_file = os.path.join(working_directory, dataset_ID + '.picked_peaks.' + mode + '.Rimage')
    pdf_file = os.path.join(working_directory, dataset_ID + '.picked_peaks.' + mode + '.pdf')
    
    ## Call R script that's a wrapper to xcsm() peak picking
    # This R file uses xcmsSet to pick peaks from the mzML files specified in mzdatafiles.<mode>.txt file, created above.
    # It saves the output diagnostics as a pdf file, and the resulting xs xcmsSet object in the rimage file.
    # It also labels each file in the xs object by its sampleID (using the sampclass method). The sample IDs are provided through the sampleIDs.<mode>.txt file, created above
    # All parameters except pdf_file and rimage_file are inputs to this function.
    cmdstr = 'Rscript /home/ubuntu/scripts/pick_peaks.R -p ' + ppm + ' -s ' + snthresh + ' --filterMin ' + prefilter_min + \
              ' --filterMax ' + prefilter_max  + ' -i ' + integrate + ' --peakMin ' + peakwidth_min + \
              ' --peakMax ' + peakwidth_max + ' -n ' + noise + ' -f ' + mzdatafiles + \
              ' --pdf ' + pdf_file + ' --rimage ' + rimage_file + ' --sids ' + sidsfile
    os.system(cmdstr)
    
    ## Track processing that's been done on samples by writing a new sequence file for just these samples
    tmp_seq_df = copy(seq_df[seq_df['ion mode'] == mode])
    tmp_seq_df['mzML files'] = files  
    tmp_seq_df['pick_peaks_rimage'] = len(files)*[rimage_file]
        
    proc_tracker = os.path.join(working_directory, dataset_ID + '.processing_tracker.' + mode + '.txt')
    tmp_seq_df.to_csv(proc_tracker, sep='\t')
    
    return rimage_file, proc_tracker

def align_peaks(rimage, batch, samples, mode, proc_file, working_directory, dataset_ID):
    """
    Aligns peaks for the samples in the given batch. 
    
    Inputs:
    =======
    - rimage:              (str) path to Rimage file containing xs xcmsSet object with picked peaks
    - batch:               (str) name of current batch being processed
    - samples:             (list) sample IDs of samples in batch. Should be as in the first column of the sequence file (index of seq_df)
    - mode:                (str) 'positive' or 'negative'
    - proc_file:           (str) path to file that is tracking the processing done so far
    - working_directory:   (str) path to working directory where files will be saved
        
    Creates:
    ========
    - <dataset_ID>.all_peaks.batch_<batch>.<mode>.csv        (file) comma-delimited file with all peaks for all samples in batch
    - <dataset_ID>.aligned_table.batch_<batch>.<mode>.csv    (file) comma-delimiated file with aligned peaks, isotopes, and adducts for all samples in batch
    
    """

    ## Define output file names
    all_peaks = os.path.join(working_directory, dataset_ID + '.all_peaks.batch_' + batch + '.' + mode + '.csv')
    aligned_table = os.path.join(working_directory, dataset_ID + '.aligned_table.batch_' + batch + '.' + mode + '.csv')

    ## Read in the proc_file
    proc_df = pd.read_csv(proc_file, sep='\t', index_col=0)
    proc_df['batch ' + batch] = [1 if s in samples else 0 for s in proc_df.index]
    
    ## Write in a temp file indicating whether each sample is in the batch or not. Contains a 'inbatch' column where each sample is yes (1) or no (0)
    # The xcmsSet object in R has each sample labeled by its sampleID. The sampleIDs are in the same order as in proc_df, so the classlist given to the R script should have the samples in that order as well.
    # The R script reads in this tmp_batch_index.txt file to get the classlist. This file indicates (with 1's and 0's) which samples to align
    tmp_file = os.path.join(working_directory, 'tmp_batch_index.txt')
    with open(tmp_file, 'w') as f:
        f.write('sample\tinbatch\n')
        for sid in proc_df.index:
            f.write(sid + '\t' + str(proc_df.loc[sid, 'batch ' + batch]) + '\n')            
    
    ## Call align_peaks.R to align the peaks and find adducts + isotopes
    # aligned_table and all_peaks are output files. Everything else is an input
    os.system('Rscript /home/ubuntu/scripts/align_peaks.R --rimage ' + rimage + ' --batch ' + tmp_file + ' --mode ' + mode + 
              ' --aligned ' + aligned_table + ' --allpeaks ' + all_peaks)
    
    ## Update the processing tracker file with the file name of the aligned_table (instead of 1's and 0's as above)
    proc_df['batch ' + batch] = [aligned_table if s in samples else 0 for s in proc_df.index]
    proc_df.to_csv(proc_file, sep='\t')
    
    ## Open the aligned feature table and replace columns with sample IDs (instead of file names)
    # The sample names given by xcms are the filename minus the .mzML extension
    # If mode is negative, these files will be <stuff>.threshold1000, otherwise they will just be <stuff>
    # <stuff> is in the "file name" column of the sequence file
    # Make a dict with {file name in xcms: sample ID}
    fnames = proc_df['file name']
    if mode == 'negative':
        fnames = [f + '.threshold1000' for f in fnames]
    fname2sid = {f:s for f, s in zip(fnames, proc_df.index)}    
    
    # Now replace the columns with the sample IDs
    aligned_df = pd.read_csv(aligned_table, sep=',')
    new_cols = list(aligned_df.columns)
    for i in range(0, len(new_cols)):
        if new_cols[i] in fname2sid:
            new_cols[i] = fname2sid[new_cols[i]]
    aligned_df.columns = new_cols
    print('[[Align peaks]] Saving aligned feature table as ' + aligned_table)
    aligned_df.to_csv(aligned_table, sep=',')


#%%# Other helpful functions for metabolomics preprocessing
def extract_batches(seq_df, mode):
    # Extract the batches of samples we want to align.
    # Batches should be specified in the sequence file in the column 'batches'
    # If one sample should be in multiple batches, the multiple batches should be comma-separated
    # Returns batches, a list of the batch names to process that have the right mode
    #         b2s, a dictionary with {batch: [samples in that batch]}. Samples are labeled by their sample ID in the sequence file (which are the indices of seq_df)

    # Need to keep samples in the same order as in the sequence files
    all_samples = seq_df.index

    ## Create dictionary of batches for each sample. s2b[sample] = [batches that sample is in]
    s2b = {key: value for key, value in zip(all_samples, seq_df['batches'])}
    batches = []
    for s in all_samples:
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
        b2s[batch] = [s for s in all_samples if batch in s2b[s]]
        
    ## Extract the batches that have all samples with the specified mode
    # If there's a batch with mixed sample modality, remove that batch too
    b2m = {}
    to_remove = []
    for batch in batches:
        b2m[batch] = [seq_df.loc[s, 'ion mode'] for s in b2s[batch]]
        mode_match = sum([i == mode for i in b2m[batch]])/float(len(b2m[batch]))   # this gets the percent of samples in each batch that match the input mode. Bc these are both ints, it's either 1 or 0.
        if mode_match != 1:
            to_remove.append(batch)
    for batch in to_remove:
        batches.remove(batch)
    
    return batches, b2s
