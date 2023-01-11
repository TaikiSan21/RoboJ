# Create list of all the PASCAL Wave files from stereo DASBRs, insert corrected UTC time, and truncate to times in the water

# You will need to set the working directory to the location of all
#   R scripts and data files on your computer
RootWD= "C:/Jay/ACOUSTIC/Buoy Recorder/PASCAL/Files for Submission/Supplemental R-code & data"
# This script requires access to a file with the original sound files
WAVdir= "S:/1650_PASCAL_2016/Recordings/DASBR/Original Sound Files"

# This Function will create a UTC date/time variable from PASCAL WAV filenames and correct UTC error in Soundtraps
MakeDateTimeFromFilename= function(filenames,DriftNum) {
  ntot= length(filenames)
  DateTimeStr= rep(NA,ntot)
  
  # split all file names strings at dots (periods)
  SplitStringsOnDots= strsplit(as.character(filenames),"[.]")
  
  # create matrix from list, filling NA for shorter elements of the list
  n.obs <- sapply(SplitStringsOnDots, length)  #length of each element of list
  seq.max <- seq_len(max(n.obs))               #create sequence from 1 to maximum length of list
  mat <- t(sapply(SplitStringsOnDots, "[", i = seq.max))  #create matrix from list
  
  # for Soundtrap, the split string is longer than 2 and date/time string is second element of matrix
  DateTimeStr[n.obs>2]= paste("20",mat[n.obs>2,2],sep='')
  
  # for SM2 and SM3, re-split the string based on underscore
  SplitStringsOnUnderscore= strsplit(mat[,1],"_")
  
  # create second matrix from new list, filling NA for shorter elements of the list
  n.obs2 <- sapply(SplitStringsOnUnderscore, length)  #length of each element of list
  seq.max2 <- seq_len(max(n.obs2))               #create sequence from 1 to maximum length of list
  mat2 <- t(sapply(SplitStringsOnUnderscore, "[", i = seq.max2))  #create matrix from list
  
  # for SM2 & SM3, the date/time strings are last two elements of matrix2
  DateTimeStr[n.obs2==3]= paste(mat2[n.obs2==3,2],mat2[n.obs2==3,3],sep='')
  DateTimeStr[n.obs2==4]= paste(mat2[n.obs2==4,3],mat2[n.obs2==4,4],sep='')
  
  # add date/time variable to original dataframe
  DateTime= as.POSIXct(DateTimeStr, format="%Y%m%d%H%M%S", tz="utc")
  
  # correct to true UTC time for ST4300 files that were downloaded with the wrong clock settings
  Add7Hrs_TF= (DriftNum<=6 |DriftNum==9 |DriftNum==12 |DriftNum==15 |DriftNum==16 | DriftNum>=26)
  Add10.8Hrs_TF= (DriftNum == 1)
  DateTime[Add7Hrs_TF]= DateTime[Add7Hrs_TF] + (7 * 60 * 60)
  DateTime[Add10.8Hrs_TF]= DateTime[Add10.8Hrs_TF] + (10.8 * 60 * 60)
  
  return(DateTime)
}

setwd("WAVdir")
filelist= list.files(pattern='*.wav',recursive=TRUE)

# Eliminate the single-channel soundtrap files
filelist2= filelist[regexpr(pattern=glob2rx("*201883689*"),text=filelist)<0]       #remove Soundtrap ST300HF-1 files
filelist2= filelist2[regexpr(pattern=glob2rx("*1678299174*"),text=filelist2)<0]       #remove Soundtrap ST300HF-Jason files
filelist2= filelist2[regexpr(pattern=glob2rx("*Stitched*"),text=filelist2)<0]       #remove old stitched-together soundtrap files
filelist2= filelist2[regexpr(pattern=glob2rx("*DASBR3*"),text=filelist2)<0]       #remove Station 14 with unuseable recordings
length(filelist2)

# create data frame with file name and station number
Station= as.numeric(substr(filelist2,9,10))                                       #read 2-digit station number
Station[is.na(Station)]= as.numeric(substr(filelist2[is.na(Station)],9,9))        #if NA, read just 1-digit station number
AllSoundFiles= data.frame(File=filelist2,Station=Station)

AllSoundFiles$BaseName= basename(as.character(AllSoundFiles$File))

DriftNum= as.numeric(AllSoundFiles$Station)
AllSoundFiles$UTC_DateTime= MakeDateTimeFromFilename(AllSoundFiles$BaseName,DriftNum)

# add recorder type to dataframe
AllSoundFiles$Recorder= as.character('ST4300')
AllSoundFiles$Recorder[DriftNum == 7]=  'SM3M'
AllSoundFiles$Recorder[DriftNum == 10]= 'SM3M'
AllSoundFiles$Recorder[DriftNum == 13]= 'SM3M'
AllSoundFiles$Recorder[DriftNum == 17]= 'SM3M'
AllSoundFiles$Recorder[DriftNum == 8]=  'SM2Bat'
AllSoundFiles$Recorder[DriftNum == 11]= 'SM2Bat'
AllSoundFiles$Recorder[DriftNum == 14]= 'SM2Bat'
AllSoundFiles$Recorder[DriftNum == 18]= 'SM2Bat'
AllSoundFiles$Recorder= factor(AllSoundFiles$Recorder)

# note whether good to use (after deployment date/time and before retrieval datetime)
# get deployment and retrieve times in Excel format and convert to R format
setwd(RootWD)
DeployRetrTimes= read.csv("DeplDetrDateTimes.csv")
DeployTime= as.POSIXct((DeployRetrTimes$UTC_DeployTime-2)*24*60*60,origin='1900-01-01',tz='utc')
RetrTime= as.POSIXct((DeployRetrTimes$UTC_RetrTime-2)*24*60*60,origin='1900-01-01',tz='utc')

# determine whether each file is good (within deploy and retrieve time)
Good2Use= rep(FALSE,length(AllSoundFiles$Recorder))
for (iter in 1:30) {
  Good2Use[DriftNum==iter]= (AllSoundFiles$UTC_DateTime[DriftNum==iter] >= DeployTime[iter]) &
                                          (AllSoundFiles$UTC_DateTime[DriftNum==iter] <= RetrTime[iter])
}
AllSoundFiles= AllSoundFiles[Good2Use,]

# output new dataframe with added variables
setwd("RootWD")
write.csv(AllSoundFiles,file="PASCAL SoundFileList.csv")

