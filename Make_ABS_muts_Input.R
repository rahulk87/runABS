## Run to make the mutation input for ABSOLUTE

## setwd('/Users/selenicp/Documents/Projects/LCIS/LCIS Exonseq/')

## You can change this to any file for where you stored the mutation summary sheet. 
mutationsFile <- "muts.csv" 
outRData <- "mutations.rda"

mutationSummary <- read.csv(mutationsFile, header = TRUE, stringsAsFactors = FALSE)
sampleList <- unique(mutationSummary$TUMOR_SAMPLE)

## create a data.frame containing the colums required to run ABSOLUTE
## columns are in the same order as the `cbind.data.frame` previously used
## Tumor_Sample_Barcode, Hugo_Symbol, t_ref_count, t_alt_count, dbSNP_Val_Status, Chromosome, Start_position
absIn <- data.frame(Tumor_Sample_Barcode = as.character(mutationSummary$TUMOR_SAMPLE),
                    Hugo_Symbol = paste(mutationSummary$ANN....GENE, mutationSummary$ANN....HGVS_P, sep = "_"),
                    t_ref_count = mutationSummary$TUMOR_DP - round(mutationSummary$TUMOR_DP * mutationSummary$TUMOR_MAF),
                    t_alt_count = round( mutationSummary$TUMOR_DP * mutationSummary$TUMOR_MAF),
                    dbSNP_Val_Status = rep("validated", nrow(mutationSummary)),
                    Chromosome = as.numeric(gsub("X", "23", as.character(mutationSummary$CHROM))),
                    Start_position = mutationSummary$POS)


## writes a .txt file per tumor sample
sapply(sampleList, function(TM) write.mutationSummaryle(absIn[absIn$Tumor_Sample_Barcode == TM , ], file = paste(TM, "muts.txt", sep = "_"), sep = "\t", row.names = FALSE))

save.image(outRData)
