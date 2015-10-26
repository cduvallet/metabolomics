#!/usr/bin/python

import preprocessing_metab as mtab

working_directory = '/home/claire/Dropbox\ \(MIT\)/MIT/alm_lab/metabolomics'
neg_files = working_directory + '/test_data/test_fnames_neg.txt'

neg_params = {}
neg_params['ppm'] = 2
neg_params['snthresh'] = 10
neg_params['prefilter_min'] = 5
neg_params['prefilter_max'] = 1000
neg_params['integrate'] = 2
neg_params['peakwidth_min'] = 20
neg_params['peakwidth_max'] = 60
neg_params['noise'] = 1000

mtab.pick_peaks(neg_files, 'neg', neg_params, working_directory)
