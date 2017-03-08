## Run to make the mutation input for ABSOLUTE
## setwd('/Users/selenicp/Documents/Projects/LCIS/LCIS Exonseq/')
filename    = 'muts.csv'     # You can change this to any file for where you stored the mutation summary sheet. 

tab         = read.csv(filename, header=T, stringsAsFactors = F)
sample_list =  unique(tab$TUMOR_SAMPLE)

## This goes through and removes the necessary columns from the mutation summary and 
## reformat it to be the correct input for ABSOLUTE
for (i in 1:length(sample_list)){
  sample = sample_list[i]
    Tumor_Sample_Barcode = as.character(tab$TUMOR_SAMPLE)
    Hugo_Symbol          = paste(tab$ANN....GENE ,tab$ANN....HGVS_P,sep="_")
    Tumor.Maf            = as.numeric(tab$TUMOR_MAF)
    x                    = tab$TUMOR_DP*Tumor.Maf
    t_alt_count          = round(x,0)
    t_ref_count          = tab$TUMOR_DP - t_alt_count

    for (i in nrow(tab)) {
        dbSNP_Val_Status = c(rep("validated",i))
    }
    
    Chromosome         = as.character(tab$CHROM)
    Chromosome         = as.numeric(replace(Chromosome, Chromosome=="X", 23)) #Change chrX to chr23
    Start_position     = as.numeric(tab$POS)
    bind<-cbind.data.frame(Tumor_Sample_Barcode,Hugo_Symbol,t_ref_count,t_alt_count,dbSNP_Val_Status,Chromosome,Start_position)
    bind = subset(bind, subset=bind$Tumor_Sample_Barcode == sample)
    write.table(bind, paste(sample, 'muts.txt', sep="_"), sep="\t",row.names = F) # write separate file for each sample
}
