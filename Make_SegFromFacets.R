# Makes ABS seg file from Facets seg file

# load funcions library
source("fxLib.R")

#setwd('/Users/selenicp/Documents/Projects/Radiogenomics/WithinPatientSufam/ABS/')  # If you run this directly after Make_ABS_muts_Input.R you dont need change the directory

cncfDir <- "~/mskcc/BioInformatics/Leiomyoma_Progression/"
files <- list.files(path = ".", pattern = "cncf.txt$") # List all segments files

## consider alterantive :
##  sapply(file.path(cncfDir, files), make_seg)

ABS_filename = as.numeric(0)

for (i in 1:length(files)){
    make_seg(file.path(cncfDir, files[i]))

  ABS_filename[i] = substr(files[i],1,regexpr("_",files[i])-1)

}
