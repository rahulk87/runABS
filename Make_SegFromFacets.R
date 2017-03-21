# Makes ABS seg file from Facets seg file
# load funcions library
source("fxLib.R")

load("mutations.rda")

cncfDir <- "cncf/"
files <- list.files(path = cncfDir, pattern = "cncf.txt$") # List all segments files

##  output files are placed on the current directory
##  to have output files in the cncf/ directory, make the files object contain the full path i.e. file.path(cncfDir, file.cncf.txt), "."
sapply(files, make_seg, cncfDir)

save.image("mutations.rda")
