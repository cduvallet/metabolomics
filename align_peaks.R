
library(optparse)
library(xcms)
library(CAMERA)

## Uses xcms to align peaks in a batch of samples
option_list = list(
  make_option(c('-i', '--rimage'), default='', type='character',
              help='path to rimage file for aligning'),
  make_option(c("-b", "--batch"), default='', type='character',
              help='path to file indicating which samples to aling. File should have sample names from rimage in one column and inbatch in the other (either 0 or 1)'),
  make_option(c('--allpeaks'), default='', type='character',
              help='path to file containing all the peaks in individual samples'),
  make_option(c('-o', '--aligned'), default='', type='character',
             help='path to file with the aligned feature table'),
  make_option(c('-m', '--mode'), default='', type='character',
              help='ionization mode. For use in finding adducts')
)


args = parse_args(OptionParser(option_list=option_list), args=commandArgs(TRUE))

## Read in which samples are part of this processing batch
#print(args$rimage)
#print(args$batch)
classes = read.delim(args$batch, stringsAsFactors=FALSE)$inbatch
#print(classes)

# For some reason, loading an Rimage erases my args (but doesn't erase these variables)
alignedfile = args$aligned
allpeaksfile = args$allpeaks
mode = args$mode
## Load Rimage which contains xs (xcmsSet object)
load(args$rimage)

## Take just the subset of xs for the samples in this batch
sampclass(xs) = classes
xsSubset = split(xs, classes)[['1']]
nsamples = length(xsSubset@filepaths)
#print(xsSubset)

## Align peaks
# Align across samples with Obiwarp. Returns an xcmsSet object. I have no idea how rc.obi and xsSubset are different from each other
# All of these hard-coded parameters are from Krista's script. Not sure if they're worth fiddling with
rc.obi = retcor.obiwarp(xsSubset, plottype='deviation', profStep=0.1, distFunc='cor', gapInit=0.3, gapExtend=0.4)
# Now group peaks from different samples together
xgN3 = group.density(xsSubset, minfrac=0, minsamp=1, bw=30, mzwid=0.001)
# And fill up any peaks that weren't found in individual samples but are considered real peaks in some other samples
xgF = fillPeaks(xgN3, method='chrom')

## I don't know why she removes all of these objects...?
## Should at some point figure out if the retcor.obiwarp and group.density functions modify xsSubset in place or not
## xgF is the final xcmsSet object with all the stuff we want.
#rm(rc.obi, xgN3, xsSubset)

## Write current allpeaks and aligned_peaks files
write.csv(xgF@peaks, file=allpeaksfile)
write.csv(peakTable(xgF), file=alignedfile)

## Now do CAMERA - code copied from Krista Longnecker
xsa = xsAnnotate(xgF)

#group the features initially just by retention time
xsaF = groupFWHM(xsa)

#figure out which features also have a matching 13C feature. Have to enter both the relative error (ppm) and the absolute error (mzabs)
xsaFI = findIsotopes(xsaF, ppm=1.5, mzabs = 0.0001, minfrac = 1/nSamples)

#now group by the correlations based on (1) intensity, (2) EIC, (3) isotopes...
xsaC = groupCorr(xsaFI, cor_eic_th=0.75, pval=0.05, graphMethod="hcs",
                 calcIso=TRUE, calcCiS=TRUE, calcCaS=FALSE)

#setup the file to also look for adducts, only go with the primary adducts for the moment
if (mode == 'negative'){
  file = system.file("rules/primary_adducts_neg.csv", package="CAMERA")
} else if (mode == 'positive'){
  file = system.file('rules/primary_adducts_post.csv', package='CAMERA')
}
  
rules = read.csv(file)
an = findAdducts(xsaC, polarity=mode, rules=rules, ppm=1.5)

#do some housecleaning
rm(xsa,xsaF,xsaFI,xsaC)

## Write final file outputs
write.csv(getPeaklist(an), file=alignedfile) 
#save.image()