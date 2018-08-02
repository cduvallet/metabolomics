Data processing code
====================

raw2feats.py
------------

``raw2feats.py`` is the main workhorse of the MS processing code. It
reads in the summary file, picks peaks (if applicable), and aligns peaks
for all batches specified in the sequence file. It basically coordinates
all of the inputs and outputs and calls wrappers to R functions as
necessary. Most of the functions that actually do work (i.e. pick and
align peaks) are found in ``preprocessing_mtab.py``.

Parsing the summary file
------------------------

The ``SummaryParserMtab.py`` module simply reads in ``summary_file.txt``
and stores its attributes in a dictionary. ``SummaryParserMtab.py``
looks for the ``summary_file.txt`` in the input directory.

Picking peaks
-------------

If an Rimage file is specified in ``summary_file.txt``, this part is
skipped. If not, ``raw2feats.py`` calls the ``pick_peaks`` function in
the ``preprocessing_mtab.py`` module. ``pick_peaks()`` calls
``pick_peaks.R`` and saves the PDF and Rimage files resulting from the
call to ``xcmsSet``. The Rimage file contains an ``xcmsSet`` object
called ``xs``. Once peaks are picked, ``summary_file.txt`` is updated
with the respective ``RIMAGE`` file so that future processing calls skip
the time-consuming peak picking step and go straight to aligning.

Aligning peaks
--------------

After peaks are picked, ``raw2feats.py`` reads in all of the specified
batches in the ``batches`` column in the sequence file. One sample may
be in multiple batches - batch names should be comma-separated in the
sample’s cell in the ``batches`` column. If a batch contains samples of
multiple ionization modes, that batch is thrown out and never processed.

``align_peaks.R`` first loads in the Rimage file that was either
specified in the summary file or created by picking peaks. It identifies
which samples to align and uses the xcms functions to align peaks across
samples, group these peaks together, and fill in any peaks that weren’t
found in individual samples but are considered real peaks in some other
samples. ``align_peaks.R`` then finds isotopes and adducts (for the
specified mode) using ``CAMERA``.

Back in ``preprocessing_mtab.py``, the sample IDs in the aligned table,
which are currently the mzML file names, are replaced by their sample ID
in the sequence file.
