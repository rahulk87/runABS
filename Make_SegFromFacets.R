# Makes ABS seg file from Facets seg file

make_seg<-function(filename){
    seg<-read.delim(filename,header=T,as.is=T)
    Chromosome    = seg$chrom
    Start         = seg$loc.start
    End           = seg$loc.end
    Num_Probes    = seg$num.mark
    Segment_Mean  = seg$cnlr.median
    seg.df        = cbind(Chromosome,Start,End, Num_Probes, Segment_Mean)
    colnames(seg.df) = c("Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
    title         = substr(filename,1,regexpr("_",filename)-1)
    write.table(seg.df,file=paste(title,".txt",sep=""),sep="\t",row.names = F)
}

#setwd('/Users/selenicp/Documents/Projects/Radiogenomics/WithinPatientSufam/ABS/')  # If you run this directly after Make_ABS_muts_Input.R you dont need change the directory
files <- list.files(pattern = "\\cncf.txt$") # List all segments files
ABS_filename = as.numeric(0)
for (i in 1:length(files)){
  make_seg(files[i])
  ABS_filename[i] = substr(files[i],1,regexpr("_",files[i])-1)
}
