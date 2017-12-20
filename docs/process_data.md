## 4 Processing your data

### 4.1 Runningraw2feats.py

To process your data, navigate to the directory containing your summary and
sequence files using the command line. If this directory does not also contain your data, make sure to specify that in the `DATA_DIRECTORY` attribute in the summary file. Then, type

```
python /path/to/raw2feats.py -i /path/to/input/directory
```

Alternatively, if you want to specify a different output directory than the default, you can type

```
python /path/to/raw2feats.py -i /path/to/input/directory -o /full/path/to/output/directory
```

Note that you must have `raw2feats.py`, `SummaryParserMtab.py`, and
`preprocessingmtab.py` all added to your `PYTHONPATH` or in the same folder.

### 4.2 Knowing when it’s finished

When the data is finished processing, you should see `[[Processing mzML files]]
Done.` print to screen.

### 4.3 Saving the outputs to files

Using Windows DOS command line, you can save the `stdout` (i.e. anything
that is spit out by any of the python codes) by piping your command to a file
name using `>`. You can also save the outputs of the calls to the R files using
`2>`.
For example, the following command saves the python outputs to `out1.txt`
and the R outputs to `out2.txt`

```
python /path/to/raw2feats.py -i /path/to/input/directory > out1.txt 2> out2.txt
```

Saving outputs on UNIX command line should be straightforward (but to be honest I wrote most of this code back in the dark ages before I got a Mac so I haven't tested it out).
