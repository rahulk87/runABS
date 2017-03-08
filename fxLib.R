## functions library

make_seg <- function(filename) {
    
    seg = read.delim(filename, header = TRUE, as.is =TRUE)

    seg.df = seg[ , c("chrom", "loc.start", "loc.end", "num.mark", "cnlr.median")]

    colnames(seg.df) = c("Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

    title = substr(filename, 1, regexpr("_", filename)-1)
    
    write.table(seg.df, file = paste(title, ".txt", sep = ""), sep = "\t", row.names = FALSE)
} # end make_seg


## unnecesary function
#ABS <- function(sampleID, primary.disease = NULL, min.ploidy = 0.95, max.ploidy = 4.5,
#                max.sigma.h = 0.7, platform = c("Illumina_WES", "SNP_6.0"),
#                copy_num_type = "total", sigma.p = 0, results.dir = "Output/",
#                max.as.seg.count = 1500, min.mut.af = 0, min_probes = 10,
#                max.non.clonal = 0, max.neg.genome = 0, ...) {
#    require(ABSOLUTE) || stop("ABSOLUTE not installed")
#    RunAbsolute(seg.dat.fn = paste(sampleID, '.txt', sep = ""),
#                min.ploidy = max.ploidy,
#                maf.fn = paste(sampleID, '_muts.txt', sep = ""), 
#                max.ploidy = max.ploidy,
#                max.sigma.h = max.sigma.h, 
#                platform = platform,
#                copy_num_type = copy_num_type,
#                sigma.p = sigma.p,
#                results.dir = results.dir,
#                primary.disease = primary.disease,
#                sample.name = sampleID,
#                max.as.seg.count = max.as.seg.count,
#                min.mut.af = min.mut.af,
#                min_probes = min_probes,
#                max.non.clonal = max.non.clonal,
#                max.neg.genome = max.neg.genome)
#} # end ABS


