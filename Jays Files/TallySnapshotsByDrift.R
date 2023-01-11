# This code will read a file with DASBR filenames from PASCAL and output file count by drift

# You will need to set the working directory to the location of all
#   R scripts and data files on your computer
RootWD= "C:/Jay/ACOUSTIC/Buoy Recorder/PASCAL/Files for Submission/Supplemental R-code & data"
AllSoundFiles= read.csv(file="PASCAL SoundFileList.csv")
AllSoundFiles$DriftNum= as.numeric(AllSoundFiles$Station)

# tally number of good files by station number
CountByDrift= rep(NA,28)
for (i in 1:28) {CountByDrift[i]= length(AllSoundFiles$DriftNum[AllSoundFiles$DriftNum==i])}
CountByDrift[14]= 0                                        #no good recordings from Station 14
cat("total file counts by drift: ",CountByDrift[1:22])
sum(CountByDrift[1:22])

outputDF= data.frame(Drift=1:22,Snapshots=CountByDrift[1:22])
write.csv(outputDF,file="SnapshotTallyByDrift.csv")
