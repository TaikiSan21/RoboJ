######## ReadPamGuardBinaries      code from Taiki Sakai
# This version is for PamGuard v2_00_16e

# Get Taiki's Packages "PamBinaries" and "PAMmisc" from Github.  
#   un-comment the following 4 lines and run to install them with the devtools package.
#   first need to install devtools and remotes from CRAN repository
# 
    # library(devtools)
    # library(remotes)
    # devtools::install_github('TaikiSan21/PamBinaries')
    # install.packages("remotes")
    # remotes::install_github("TaikiSan21/PAMmisc")

 
library(PamBinaries)    #package 'PamBinaries' requires R >=  3.4.0
library(PAMmisc)    #package 'PamBinaries' requires R >=  3.4.0
library(dplyr)
timestamp()

# Change the folder name to match the location of the binaries for a single drift or a single SM3M SD card
# folder= "P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-7\\Card-A"
# folder= "P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-7\\Card-C"
# folder= "P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-7\\Card-D"
# folder= "P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-7\\Card-AC"
# folder= "P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-17\\Card-M"
# folder= "P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-17\\Card-N"
# folder= "P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-17\\Card-O"
# folder= "P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-17\\Card-P"
# Read base names for templates

TemplNames=c("Zc","BW43","BW39V","Ms","Bb","BW70")
numTemplates= length(TemplNames)


 folder= "P:\\PASCAL\\Binaries&DBs by Drift ReRun w 2_00_16e\\Drift-22"
 numTemplates=0
 
# 
setwd(folder)


# 
# Get the names for all the click detector binaries file in the given folder and subfolders
  files= list.files(folder,pattern=glob2rx("Click_Detector_Click_Detector_Clicks*.pgdf"),recursive=TRUE,full.names=TRUE)
  nfiles= length(files)

#
# Get binary data for all files
  if (nfiles == 0) {
    cat(" ERROR, no binary click files found  \n")
    stop()
  } else {
    DFlist<- vector('list', length = nfiles)  #create list of dataframes for each binary click file
    # cycle through all binary click file names to extract binary data
    for (ifile in 1:nfiles) {
      cat(ifile," of ",nfiles,"  ",files[ifile],"\n")
      # loadPamguardBinaryFile will load any binary file type
      binaryData <- loadPamguardBinaryFile(files[ifile],skipLarge=TRUE)
      nClicks <- length(binaryData$data)
      if (nClicks > 0) {
        #create click dataframe from binary data including the click template "delta" value
        if (numTemplates==0) {
          clickDF<- pbToDf(binaryData) 
        } else {
          clickDF<- pbToDf(binaryData,TemplNames) 
        }
        clickDF$UID= as.character(clickDF$UID)
        clickDF$BinaryFile= files[ifile] #add element to click dataframe with full file name 
        DFlist[[ifile]]= clickDF      #append AnnotationDF columns to clickDF
      }
    }
  }

  cat(" Begin binding list of dataframes \n")
  MasterClickDF= bind_rows(DFlist)
  write.csv(MasterClickDF,"ClickBinaries.csv")
  timestamp()

