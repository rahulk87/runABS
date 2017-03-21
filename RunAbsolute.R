#ABSOLUTE PIPELINE
#### In R #####
## setwd('/Users/selenicp/Documents/Projects/Radiogenomics/WithinPatientSufam/ABS/Output/')

obj.name = "new" ##  Change this to describe the tumor type
load("mutations.rda")  ## from MakeABS_Muts

## libraries
library(numDeriv)
library(ABSOLUTE)
source("fxLib.R")

## STEP 1 : run absoloute per sample
##  sampleList = unique(list$Sample.ID)

sapply(sampleList, function(SM) {
    RunAbsolute(paste(SM, ".txt", sep = ""),
                min.ploidy = 0.95,
                max.ploidy = 4.5,          ## number based on copy number plot
                primary.disease = obj.name,
                max.sigma.h = 0.7,
                ## platform if you are using data from Facets use "Illumina_WES" or "SNP_6.0" from oncoscan array
                platform = "Illumina_WES",
                maf.fn = paste(SM, '_muts.txt', sep = ""), 
                copy_num_type = "total",   ## if you are using segmentation files (both from sequencing or array)
                sigma.p = 0, 
                results.dir = "Output/",   ## output path
                sample.name = SM,
                max.as.seg.count = 1500,   ## this is the total number of segments in your segmentation file.
                ## You have to increase this number if the size of your segmentation
                ## file is greater than 1500
                min.mut.af = 0,            ## this is the minimum mutation allele frequency.
                ## Mutations with lower allelic fractions will be filtered out before analysis
                min_probes = 10,
                max.non.clonal = 0,
                max.neg.genome = 0)
})

                                            
## STEP 2
## Go to the output of STEP 1 
setwd("Output/")
absolute.files = list.files(pattern = "*.RData")
indv.results.dir = "absolute_results" 
copy_num_type = "total"

CreateReviewObject(obj.name = obj.name, 
                   absolute.files = absolute.files, 
                   indv.results.dir = indv.results.dir, 
                   copy_num_type = copy_num_type, 
                   plot.modes = TRUE, 
                   verbose = TRUE)

##  Open the .PP-calls_tab.txt and pick the corect solution for each 
## STEP 3
calls.path = file.path("absolute_results", paste(obj.name, ".PP-calls_tab.txt", sep=""))
modes.path = file.path("absolute_results", paste(obj.name, ".PP-modes.data.RData", sep=""))
output.path = file.path("absolute_results", "Output")

ExtractReviewedResults(reviewed.pp.calls.fn = calls.path, analyst.id = "test", modes.fn = modes.path, out.dir.base = output.path, obj.name = "absolute", copy_num_type = "total")
