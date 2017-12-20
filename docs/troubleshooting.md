# 5 Common confusing errors

**Error in `phenoDataFromPaths(files)`: Directory tree must be level**

This error means that there is probably a mismatch between the files specified
in your sequence file and those present in your data folder. The `pick_peaks`
code goes through all files it finds in the sequence file’s `File Name` column and tries to pick peaks for them.

**A value is trying to be set on a copy of a slice from a DataFrame**

This is a non-fatal error. Don’t worry about it.

**A subdirectory or file already exists, or other copying files errors**

Errors related to making directories and copying files should also be non-fatal.
As long as the processing continues after them, you can ignore them.

**IOError: [Errno 13] Permission denied: file_name**

This error means that one of the files that the code is trying to edit, like the `processing_tracker` or `summary_file` are open in another program. Close them and try again.

**Rscript, pickpeaks.R, or alignpeaks.R not found**

Make sure that the paths to the respective files are correct. These are in the `preprocessingmtab.py` module, in either (or both, in the case of the call to Rscript) `pickpeaks()` or `alignpeaks()`. If there are spaces in the file path and you are using a DOS command prompt, enclose the string that has spaces in quotation marks (").
