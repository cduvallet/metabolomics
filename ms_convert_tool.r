################ ms convert tool #############
#use this to have an r script to convert the RAW files into whatever format I want
# KL 3/27/2014 from code I found online
#this version should work its way through the two folders with the data
#modify for AAIW experiment during DeepDOM, KL 4/10/2014

##remember that the pattern in the file name will need to be changed in addition to the file location.

here<- getwd()
 
#note that for this to work, msconvert must be in a folder that has no spaces in the file name
msconvert <- c("C:/MSConvert/msconvert.exe")

folders <- c("C:/Users/Claire/Dropbox (MIT)/MIT/alm_lab/UW/30_and_52_again")


for (ii in 1:length(folders)) {
    
setwd(folders[ii])
FILES <- list.files(recursive=FALSE, full.names=TRUE, pattern="\\_alm_")

#Notes on the 'filters' in msconvert:
#this filter must be first; this is the peak picking (convert to centroid):--filter \"peakPicking true 1-\"
#this filter will only keep peaks > 1000 (absolute intensity) --filter \"threshold absolute 1000 most-intense\"
  
  for (i in 1:length(FILES)) {
    #this is what needs to happen for negative ion mode data...allows better peak picking with
    #threshold set in MSconvert
    system(paste(msconvert, "--mzML --filter \"peakPicking true 1-\" --filter \"threshold absolute 1000 most-intense\" -o mzML_threshold1000", FILES[i]))
    
    #this will allow peak picking for MSn levels 1 and higher...needed for xcmsFragments. No threshold here 
    #because then you will lose some of the fragments
    system(paste(msconvert, "--mzML --filter \"peakPicking true 1-3\" --filter \"msLevel 1-3\" -o mzML_noThreshold", FILES[i]))
    
  }
  
  rm(i,FILES)

}
rm(ii)
setwd(here)

#do some housecleaning
rm(here,msconvert)

# i<-32
# system(paste(msconvert, "--mzXML --filter \"peakPicking true 1-\" --filter \"threshold absolute 1000 most-intense\" -o mzML", FILES[i]))
