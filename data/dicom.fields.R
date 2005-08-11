## DICOM_VRfields.pm ver 0.3
## Andrew Crabb (ahc@jhu.edu), May 2002.
dicom.VRfields <- read.table("dicom.VRfields.txt", header=TRUE,
                             colClasses=c("character","character",
                               "numeric","numeric"))

dicom.fields <- read.table("dicom.fields.txt", header=TRUE,
                           colClasses=c("character"))

