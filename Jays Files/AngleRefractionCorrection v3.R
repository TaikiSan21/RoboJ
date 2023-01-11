# ApplyRefractionCorrectionToAllEventAngles.r
# Corrections are based on a look-up table from a MATLAB raytrace program

# You will need to set the working directory to the location of all
#   R scripts and data files on your computer
setwd("C:/Jay/ACOUSTIC/Buoy Recorder/PASCAL/Files for Submission/Supplemental R-code & data")

inputFile=  "AllEventClicks.csv"
outputFile= "AllEventClicks wRC.csv"
AllEventClicksWithAngle= read.csv(file=inputFile)

RefractionCorrections= read.csv("RefractedAngleCorrections.csv",header=FALSE)
names(RefractionCorrections)= (paste("Drift-",as.character(1:30)))
Angles= 5:85  #declination angles (from horizontal) must match what was used in RayTrace program to generate 
nAngles= length(Angles)

RefractionCorrections= Angles-RefractionCorrections

# create indices for addressing the elements of the RefractionCorrection Matrix
DriftIndex= as.integer(AllEventClicksWithAngle$Station)

# convert original angle (re: up) to declination angle (re: horizontal) and create index 
#   for addressing RefractionCorrection Matrix
AllEventClicksWithAngle$Angle= AllEventClicksWithAngle$Angle - 90
hist(AllEventClicksWithAngle$Angle,xlab="Declination Angle (deg)")
AngleIndex= round(AllEventClicksWithAngle$Angle + 1 - Angles[1])

# if index does not fall between modeled range, assign zero correction value
RefractionCorrections[nAngles+1,]=0
AngleIndex[is.na(AngleIndex)]= nAngles+1
AngleIndex[AngleIndex<1]= nAngles+1
AngleIndex[AngleIndex>nAngles] = nAngles+1

# look up angle correction for each event angle
AngleCorrections= rep(NA,length(AngleIndex))
for (i in 1:length(AngleIndex)) {
  AngleCorrections[i]= RefractionCorrections[AngleIndex[i],DriftIndex[i]]
}

# create new  angle variable with refraction corrections between given angles
AllEventClicksWithAngle$AngleRefractCorrected= AllEventClicksWithAngle$Angle - AngleCorrections

# convert declination angle (re: horizontal) to detection angle (re: down)
AllEventClicksWithAngle$AngleRefractCorrected= 90 - AllEventClicksWithAngle$AngleRefractCorrected
AllEventClicksWithAngle$AngleRefractCorrected[AllEventClicksWithAngle$AngleRefractCorrected<0]= 0

# plot original angles vs corrected angles
plot(90-AllEventClicksWithAngle$Angle,AllEventClicksWithAngle$AngleRefractCorrected,xlab="Original Detection Angle (deg)",
     ylab="Corrected Angle (deg)",main="Angle Correction",pch=".",ylim=c(0,90),xlim=c(0,90))


# output file with new field for refrection corrected angles
write.csv(AllEventClicksWithAngle,file=outputFile)

summary(AllEventClicksWithAngle)

