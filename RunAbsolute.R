#ABSOLUTE PIPELINE
#### In R #####
## setwd('/Users/selenicp/Documents/Projects/Radiogenomics/WithinPatientSufam/ABS/Output/')
obj.name = "lm" ##  Change this to describe the tumor type

library(numDeriv)
library(ABSOLUTE)

## STEP 1
sigma.p <- 0
max.sigma.h <- 0.07
min.ploidy <- 0.95 
max.ploidy <- 4.5  ## (choose this based on copy number plot)
platform <- "Illumina_WES" ## (if you are using data from Facets)
 ##  platform <- "SNP_6.0" ## (if you are using data from oncoscan (copy number array))
max.as.seg.count <- 1500 ## (this is the total number of segments in your segmentation file. You have to increase this number if the size of your segmentation file is greater than 1500)
copy_num_type <- "total" ## (if you are using segmentation files (both from sequencing or array))
max.neg.genome <- 0
max.non.clonal <- 0
min.mut.af <- 0  ## (this is the minimum mutation allele frequency. Mutations with lower allelic fractions will be filtered out before analysis )
results.dir <- "Output/" ## (output path)
min_probes <- 10

ABS <- function(n, primary.disease) {
RunAbsolute(seg.dat.fn = paste( n, '.txt', sep = ""), 
 min.ploidy = min.ploidy, 
 maf.fn = paste( n, '_muts.txt', sep = ""), 
 max.ploidy = max.ploidy, 
 max.sigma.h = max.sigma.h, 
 platform = platform, 
 copy_num_type = copy_num_type, 
 sigma.p = sigma.p, 
 results.dir = results.dir,
 primary.disease = primary.disease,
 sample.name = n,
 max.as.seg.count = max.as.seg.count, 
 min.mut.af = min.mut.af,
 min_probes = min_probes,
 max.non.clonal = max.non.clonal,
 max.neg.genome = max.neg.genome)
}


##  sample_list = unique(list$Sample.ID)
invisible(
for (i in 1:length(sample_list)){
 ABS(sample_list[i], obj.name)
}
)

## STEP 2
###Go where the output of the the first step is (you can do this using the "setwd("your output path")" command in R) ###

setwd(results.dir)

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
# STEP 3
calls.path = file.path("absolute_results", paste(obj.name, ".PP-calls_tab.txt", sep=""))
modes.path = file.path("absolute_results", paste(obj.name, ".PP-modes.data.RData", sep=""))
output.path = file.path("absolute_results", "Output")

ExtractReviewedResults(calls.path,"test",modes.path,output.path,"absolute","total")
