# metabolomics

This repo contains scripts to process metabolomics data.

## MS Convert

- `ms_convert_tool.r` - Krista's R script that she used for converting files
- `msconvert_wrapper.py` - Claire's python script to convert files

## Processing data

- `raw2feats.py` - user-interfacing script that coordinates the processing. This is the script you should run to process your data.
- `preprocessing_mtab.py` - module with functions to process metabolomics data (most of these functions are calls to the R scripts)
- `align_peaks.R`, `pick_peaks.R` - R scripts which are called and actually do the processing
- `SummaryParserMtab.py` - module to parse summary file

- `mass-spec-processing-documentation.pdf` - documentation explaining how to use `raw2feats.py`

#### To do
- documentation on ReadTheDocs
- clean up unnecessary files
