
library(optparse)
library(xcms)

## Uses xcms to align peaks in a batch of samples
option_list = list(
  make_option(c('-i', '--rimage'), default='', type='character',
              help='path to rimage file for aligning'),
  make_option(c("-b", "--batch"), default='', type='character',
              help='path to file indicating which samples to aling. File should have sample names from rimage in one column and inbatch in the other (either 0 or 1)'),
  make_option(c('--allpeaks'), default='', type='character',
              help='path to file containing all the peaks in individual samples'),
  make_option(c('-o', '--aligned'), default='', type='character',
             help='path to file with the aligned feature table')
)


args = parse_args(OptionParser(option_list=option_list), args=commandArgs(TRUE))

## Read in which samples are part of this processing batch
#print(args$rimage)
#print(args$batch)
classes = read.delim(args$batch, stringsAsFactors=FALSE)$inbatch
#print(classes)

alignedfile = args$aligned
allpeaksfile = args$allpeaks
## Load Rimage which contains xs (xcmsSet object)
load(args$rimage)

## Take just the subset of xs for the samples in this batch
sampclass(xs) = classes
xsSubset = split(xs, classes)[['1']]
#print(xsSubset)

## Align peaks
# Align across samples with Obiwarp. Returns an xcmsSet object. I have no idea how rc.obi and xsSubset are different from each other
# All of these hard-coded parameters are from Krista's script. Not sure if they're worth fiddling with
rc.obi = retcor.obiwarp(xsSubset, plottype='deviation', profStep=0.1, distFunc='cor', gapInit=0.3, gapExtend=0.4)
# Now group peaks from different samples together
xgN3 = group.density(xsSubset, minfrac=0, minsamp=1, bw=30, mzwid=0.001)
# And fill up any peaks that weren't found in individual samples but are considered real peaks in some other samples
xgF = fillPeaks(xgN3, method='chrom')

# xgF is the final xcmsSet object with all the stuff we want.
rm(rc.obi, xgN3, xsSubset)

## Write some stuff to files
write.csv(xgF@peaks, file=allpeaksfile)
write.csv(peakTable(xgF), file=alignedfile)