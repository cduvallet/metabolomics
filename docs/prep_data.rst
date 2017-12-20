Preparing your dataset for processing
=====================================

You should create a directory where you put all of the data and
associated files described below.

Raw data
--------

The pipeline does not currently support raw data, and begins instead
with open-sourced mzMLfiles.

mzML files
----------

Currently, you need to convert your raw data files into mzML files
manually. Use MSConvert for this. Each file should be run through
MSConvert twice: once with a threshold of 1000, (named with the suffix
``*.threshold1000.mzML``), and once with no threshold (named
``*.mzML``). Negative mode peak-picking works better with the
thresholded files, whereas positive-mode and MS2 processing don’t need
any of the pre-thresholding by MSConvert.

Note that the associated `github
repository <https://github.com/cduvallet/metabolomics>`__ contains a
helper script, ``msconvert_wrapper.py`` which you can use to easily
convert raw data into mzML files. Note that some raw data formats (e.g.
Thermo Fisher) can only be converted with the Windows version of
MSConvert. If this is your case, we recommend installing
`Cygwin <http://www.cygwin.com/>`__ and using ``msconvert_wrapper.py``
from that command line.

Sequence file
-------------

The sequence file should be provided to you by whoever ran the samples
on the machine. The sequence file is essentially a mapping file
containing the names of each data file and its corresponding metadata.
Some metadata that is often present includes: instrument method, file
path, injection volume, etc.

The required parts of the sequence file are as follows:

-  **Sample ID**: The first column in your sequence file should be the
   sample IDs. This is how each sample will be labeled in all downstream
   processing.
-  **File Name**: This column should contain the raw data file name
   (without any extensions). The processing code assumes that the mzML
   files are created from these file names. For example, if you have
   ``mtabalmsample1`` in this column, the code assumes that the
   corresponding mzML files are m\ ``tabalmsample1.threshold1000.mzML``
   and ``mtabalmsample1.mzML``.
-  **Ion Mode**: This column contains the ion mode used for each sample.
   Accepted values are ``negative`` and ``positive``.
-  **Batches**: This column specifies the "batches" of samples you want
   to align. When aligning the picked features to created an aligned
   feature table, you may only want to consider a subset of your samples
   (e.g. all PPL samples in one batch, all direct injection samples in
   another). Each sample can be in many (or no) batches. If a sample is
   in multiple batches, the batch names should be comma-separated within
   the same cell. The order of batches in a cell doesn’t matter, but the
   case does. If a batch contains samples from multiple ionization
   modes, that batch will not be aligned. Each batch yields an aligned
   feature table with only the samples in that batch aligned.

The first row of the sequence file should contain at least the following
case- insensitive column headers: ``SampleID``, ``File Name``,
``Ion Mode``, and ``batches``. If there is an additional line at the top
of the file (above the column headers), delete it before providing to
the pipeline. The sequence file is assumed to be comma-separated, but
the delimitation can be specified in the summary file with the attribute
``SEQUENCE_FILE_SEPARATOR``. To specify a different sequence file
delimiter, include the Pythonic string representation of the separator.
For example, a tab-delimited sequence file would have a
``SEQUENCE_FILE_SEPARATOR`` of ``\t``.

Summary File
~~~~~~~~~~~~

Once you have your data and sequence file all sorted, you need to create
a sum- mary file that “talks” to ``raw2feats.py`` (through the
``SummaryParserMtab.py`` module). Your summary file should be a
tab-delimited file named ``summaryfile.txt`` and placed in the same
directory as your sequence file.

Required attributes
^^^^^^^^^^^^^^^^^^^

The following attributes are required to be specified in
``summaryfile.txt``:

+-------------------+---------------------------------------------------------+
| ``DATASET\_ID``   | Identifier for this processing run. Output files will   |
|                   |  contain this string as an identifier.                  |
|-------------------+---------------------------------------------------------+
| ``MODE``          | Ionization mode. Accepted values are ``negative`` or    |
|                   | ``positive``. If you have both positive and negative    |
|                   | mode files to process, you will need to do them in two  |
|                   | separate runs.                                          |
|-------------------+---------------------------------------------------------+
| ``SEQUENCE\_FILE``| Name of sequence file. Full path is not necessary, as   |
|                   | sequence file is assumed to (and should) be in the same |
|                   | directory as the summary file.                          |
+-------------------+---------------------------------------------------------+


Optional attributes
^^^^^^^^^^^^^^^^^^^

The following attributes can be specified in the summary file, but are
not required for processing:

+---------------------------------+--------------------------------------------+
| ``DATA\_DIRECTORY``             | If the data is in a different directory,   |
|                                 | you can provide the full path to the       |
|                                 | directory here. Otherwise the code         |
|                                 | assumes that all of your mzML files are in |
|                                 | the input directory.                       |
+---------------------------------+--------------------------------------------+
| ``SEQUENCE\_FILE\_DELIMITER``   | The sequence file is assumed to be         |
|                                 | comma-delimited. If this is not the case,  |
|                                 | specify the delimitation here. (i.e. if    |
|                                 | your sequence file is tab-delimited, this  |
|                                 | should be ``\t``).                         |
+---------------------------------+--------------------------------------------+
| ``RIMAGE``                      | If peak-picking has already been performed |
|                                 | on this dataset, you may provide the full  |
|                                 | path to an Rimage file containing this     |
|                                 | results to load up and skip the            |
|                                 | peak-picking. This Rimage should contain   |
|                                 | an ``xcmsSet`` object named ``xs``. If no  |
|                                 | Rimage file is specified, this attribute   |
|                                 | will be updated with the correct file once |
|                                 | peak-picking has been run once, so you may |
|                                 | re-use this summary file to run different  |
|                                 | alignments without necessarily re-picking  |
|                                 | peaks.                                     |
+---------------------------------+--------------------------------------------+
| ``RAW\_DATA``                   | ``True`` if you are providing raw data     |
|                                 | that needs to be converted to mzML,        |
|                                 | ``False`` if you are directly providing    |
|                                 | mzML files. **Note: the current code does  |
|                                 | NOT accept raw data.** This functionality  |
|                                 | may be added in future iterations.         |
+---------------------------------+--------------------------------------------+

Metabolomics summary file attributes should be located between
``#mtab_start`` and ``#mtab_end`` in the summary file.

Sample Summary File
^^^^^^^^^^^^^^^^^^^

Note: all blank spaces are tab characters.

::

    DATASET_ID test

    #mtab_start
    MODE negative
    SEQUENCE_FILE test_sequence_file.csv
    #mtab_end

Sample Summary File, with optional attributes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    DATASET_ID test

    #mtab_start
    MODE negative
    SEQUENCE_FILE testsequencefile.csv
    RIMAGE full/path/to/file.Rimage
    SEQUENCE_FILE_DELIMITER \t
    #mtab_end
