
library(optparse)
library(xcms)

# Uses the xcmsSet function to pick peaks from all files given
# Takes as input the following parameters

option_list = list(
  make_option(c("-p", "--ppm"), default=2, type='integer',
              help='ppm error to use in xcmsSet'),
  make_option(c('-s', '--snthresh'), default=10, type='integer',
              help='signal to noise threshold for use in xcmsSet'),
  make_option(c('--filterMin'), default=5, type='integer',
              help='pre-filter min value'),
  make_option(c('--filterMax'), default=1000, type='integer',
              help='pre-filter max value'),
  make_option(c('-i', '--integrate'), default=2, type='integer',
              help='value for integrate parameter in xcmsSet'),
  make_option(c('--peakMin'), default=20, type='integer',
              help='peakWidth min value'),
  make_option(c('--peakMax'), default=60, type='integer',
              help='peakWidth max value'),
  make_option(c('-n', '--noise'), default=1000, type='integer',
              help='noise value for xcmsSet'),
  make_option(c('-f', '--fnames'), default='', type='character',
              help='path to file containing full paths to files to be run through xcsmSet'),
  make_option(c('--pdf'), default='picked_peaks.pdf', type='character',
              help='name of pdf file to write CentWave peak picking diagnostics to'),
  make_option(c('--rimage'), default='picked_peaks.Rimage', type='character',
                help='name of Rimage file to save picked peaks')
)

args = parse_args(OptionParser(option_list=option_list), args=commandArgs(TRUE))

### Need full path to these data files...
mzdatafiles <- read.csv(args$fnames, stringsAsFactors=FALSE)$fnames

## How many CPU cores has your machine (or cluster) ?
nSlaves=4

pdf(file = args$pdf) #if want PDF file, add this after mzdiff: ,sleep = 0.0001 and uncomment dev.off()

#For negative ion mode: ppm = 2 seems best
#For positive ion mode: ppm = 3 seems best
xs<-xcmsSet(mzdatafiles, method="centWave", ppm=args$ppm ,snthresh=args$snthresh,
            prefilter =  c(args$filterMin, args$filterMax), mzCenterFun="wMean",integrate=args$integrate, 
            verbose.columns=TRUE, peakwidth=c(args$peakMin,args$peakMax), fitgauss=TRUE, noise=args$noise, 
            mzdiff=-0.005, nSlaves=nSlaves, sleep=0.00001)
dev.off()

# Define the samples in xs$phenoData using sampclass(xs)
sampclass(xs) = mzdatafiles

save.image(args$rimage)
