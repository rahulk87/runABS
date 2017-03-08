## Run to make the mutation input for ABSOLUTE

## setwd("./")

## You can change this to any file for where you stored the mutation summary sheet
mutationsFile <- "muts.csv" 
outRData <- "mutations.rda"

## read in the mutation summary file, and create a list with the unique samples.
## fields containing NA, a period '.' or empty are treaded as missing values
mutationSummary <- read.csv(mutationsFile, header = TRUE, stringsAsFactors = FALSE, na.strings = c(NA, ".", ""))
sampleList <- unique(mutationSummary$TUMOR_SAMPLE)

## create an empty data.frame containing the number of colums required to run ABSOLUTE
## columns are in the same order as the `cbind.data.frame` previously used
## Tumor_Sample_Barcode, Hugo_Symbol, t_ref_count, t_alt_count, dbSNP_Val_Status, Chromosome, Start_position
absIn <- data.frame(matrix(NA, nrow = nrow(mutationSummary), ncol = 7))
names(absIn) <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "t_ref_count", "t_alt_count", "dbSNP_Val_Status", "Chromosome", "Start_position")

## populate absIn data.frame
absIn$Tumor_Sample_Barcode <- as.character(mutationSummary$TUMOR_SAMPLE)
absIn$Hugo_Symbol <- paste(mutationSummary$ANN....GENE, mutationSummary$ANN....HGVS_P, sep = "_")
absIn$t_alt_count <- round( mutationSummary$TUMOR_DP * mutationSummary$TUMOR_MAF)
absIn$t_ref_count <- mutationSummary$TUMOR_DP - absIn$t_alt_count
absIn$dbSNP_Val_Status <- rep("validated", nrow(mutationSummary))
absIn$Chromosome <- as.numeric(gsub("X", "23", as.character(mutationSummary$CHROM)))
absIn$Start_position <- mutationSummary$POS


## loop trhough sampleList and write a .txt file per sample
sapply(sampleList, function(TM) write.table(absIn[absIn$Tumor_Sample_Barcode == TM , ], file = paste(TM, "muts.txt", sep = "_"), sep = "\t", row.names = FALSE, na = "."))

save.image(outRData)
