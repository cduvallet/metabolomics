.. untargeted-metabolomics-pipeline documentation master file, created by
   sphinx-quickstart on Tue Dec 19 22:00:04 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Untargeted metabolomics processing pipeline
============================================================

This documentation describes the process of going from open-source mzML mass
spectrometry files to aligned feature table(s). This is orchestrated by the script `raw2feats.py`. Each dataset folder must contain a machine-readable text file called a summary file which gives instructions to `raw2feats.py`. The format of this summary file, and all files required for MS processing, are described in this documentation.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   prep_data
   proc_code
   process_data
   troubleshooting



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
