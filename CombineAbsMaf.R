## 

library(qdap)

MAF.directory <- "absolute_results/Output/reviewed/SEG_MAF/"
## results.directory <- "../../../../"

options(scipen = 999) # Sets options to not write anything in Scientific notation (This messes up Clonal Status)


## set working directory
setwd(results.directory)

total.muts <- read.csv("muts.csv", header = TRUE, stringsAsFactors = FALSE)
 # total.muts <- tab
total.muts$CHROM <- as.numeric(replace(total.muts$CHROM, total.muts$CHROM == "X", 23))
x <- as.numeric(as.character(total.muts$TUMOR_DP))*total.muts$TUMOR_MAF
total.muts$t_alt_count <- round(x, 0)
total.muts$t_ref_count <- as.numeric(as.character(total.muts$TUMOR_DP)) - total.muts$t_alt_count


## Add Ref and Alt counts to mutation file
setwd(MAF.directory)

( file_list <- list.files(pattern = "ABS_MAF.txt$") )
i <- 1
rm(ABS.MAF)


# ABS.MAF <- do.call(rbind, lapply(file_list, read.table, header = TRUE, sep = "\t"))

for (file in file_list) {
    
    if (!exists("ABS.MAF")) { # if the merged ABS.MAF doesn't exist, this creates it

        ABS.MAF <- read.table(file.path(MAF.directory, file), header = TRUE, sep = "\t")
        ABS.MAF$sample = file

    } else {

        temp_ABS.MAF <- read.table(file, header = TRUE, sep = "\t")
        temp_ABS.MAF$sample <- file
        ABS.MAF <- rbind(ABS.MAF, temp_ABS.MAF)
        rm(temp_ABS.MAF)
    }

    print(i)
    i = i+1
}

## Merge CCF values with muts.csv file
## Create empty vectors
rownames(total.muts) <- c(1:nrow(total.muts))
alt2 <- vector(mode = "numeric", length = 0)
ref2 <- vector(mode = "numeric", length = 0)
Pr_somatic_clonal <- vector(mode = "numeric", length = 0)
Cancer_Cell_Fraction <- vector(mode = "numeric", length = 0)
CI95_low <- vector(mode = "numeric", length = 0)
CI95_high <- vector(mode = "numeric", length = 0)
Clonal_Status <- vector(mode = "numeric", length = 0)

## Pull necessary columns out of MAF files
for (i in 1:nrow(total.muts)) {
    comb = subset(ABS.MAF, ABS.MAF$ref == total.muts$t_ref_count[i] & ABS.MAF$Start_position == total.muts$POS[i] & ABS.MAF$alt == total.muts$t_alt_count[i])
 comb = comb[!duplicated(comb), ]
 print(paste(i,nrow(comb), sep = "_"))

 if (nrow(comb) > 0) {
     alt2[i] = comb$alt
     ref2[i] = comb$ref
     Pr_somatic_clonal[i] = as.numeric(round(comb$Pr_somatic_clonal,5))
     Cancer_Cell_Fraction[i] = comb$cancer_cell_frac
     CI95_low[i] = ifelse(!is.na(comb$ccf_CI95_low), comb$ccf_CI95_low, -1)
     CI95_high[i] = ifelse(!is.na(comb$ccf_CI95_high), comb$ccf_CI95_high, -1)

     if (Pr_somatic_clonal[i] >= 0.5 | CI95_low[i] >= 0.9) {
         Clonal_Status[i] = "Clonal"
     } else {
         Clonal_Status[i] = "Subclonal"}
     
 } else {
     alt2[i] = "."
     ref2[i] = "."
     Pr_somatic_clonal[i] = "."
     Cancer_Cell_Fraction[i] = "."
     CI95_low[i] = "."
     CI95_high[i] = "."
     Clonal_Status[i] = "."
 }
}

## Bind new MAF columns to old mutation df and Check that tables are lines up correctly
ABS_final <- cbind.data.frame(total.muts, Cancer_Cell_Fraction, Pr_somatic_clonal, CI95_high, CI95_low, Clonal_Status, alt2, ref2)
ABS_final <- subset(ABS_final, subset = ABS_final$alt2 != ".")


## Fixing AA changes
remove_pipe <- ABS_final
temp <- data.frame(matrix(ncol = ncol(remove_pipe), nrow = 1))

repeat {
    for (k in 1:nrow(remove_pipe)) {
        if (substring(remove_pipe$ANN....HGVS_P[k], 1, 2) == ".|") {
            print(k)
            
            remove_pipe$ANN....HGVS_P[k] = substring(remove_pipe$ANN....HGVS_P[k], 3, nchar(remove_pipe$ANN....HGVS_P[k]))
            remove_pipe$ANN....EFFECT[k] = char2end(remove_pipe$ANN....EFFECT[k], "|")
            remove_pipe$ANN....GENE[k] = char2end(remove_pipe$ANN....GENE[k], "|")
        }
    }
    
    temp <- subset(remove_pipe, subset = substring(remove_pipe$ANN....HGVS_P, 1, 2) == ".|")
    if (nrow(temp) < 1) {
        break
    }
}


## Redo LOH calls
## This part only works if you set the Facets cncf.txt files in the directory

setwd("../../../../../")

muts <- remove_pipe

Positions <- cbind.data.frame(muts$TUMOR_SAMPLE, as.numeric(as.character(muts$CHROM)), as.numeric(muts$POS))
colnames(Positions) <- c("Sample.ID", "Chromosome", "Position")
Positions$Chromosome <- as.numeric(replace(Positions$Chromosome, Positions$Chromosome == "X", 23))

samp <- as.numeric(0)
samp.upd <- as.numeric(0)
Pos <- as.numeric(0)
Chr <- as.numeric(0)
tcn <- as.numeric(0)
lcn <- as.numeric(0)
loh <- as.numeric(0)

## This part is not perfect - it spits out an error if the mutation is not within one of the segments
## If it fails you go into where 100000 to the number that it gives you if it fails
## You will have to go through those by hand using the segments file
## results.directory <- "/Users/selenicp/Documents/Projects/Radiogenomics/WithinPatientSufam/ABS/"

files <- list.files(pattern = "cncf.txt$")

for (i in 1:nrow(Positions)) {
    print(i)
    samp[i] = as.character(Positions$Sample.ID[i])
    FACET.df = read.table(files[as.data.frame(strsplit(files, "_"))[1,] == samp[i]], header = TRUE)
    median.logr = median(FACET.df$mafR)
    std.logr = sd(FACET.df$mafR)
    ## FACET.df = read.table(files[substr(files, 1, 12) == samp[i]], header = TRUE)
    Chr[i] = Positions$Chromosome[i]
    Pos[i] = Positions$Position[i]
    if (i %in% c(310, 362)) {
        tcn[i] = "error"
        lcn[i] = "error"
        loh[i] = "error"
    } else {
        tcn[i] = as.character(FACET.df$tcn.em[FACET.df$chrom == Chr[i] & FACET.df$loc.start<Pos[i] & FACET.df$loc.end>Pos[i]])
        lcn[i] = as.character(FACET.df$lcn.em[FACET.df$chrom == Chr[i] & FACET.df$loc.start<Pos[i] & FACET.df$loc.end>Pos[i]])
        loh[i] = ifelse(tcn[i]>0, ifelse(lcn[i]>0, ".", "loh"), "error")
    }
    if(is.na(loh[i])){
        loh[i] = ifelse(FACET.df$mafR[FACET.df$chrom == Chr[i] & FACET.df$loc.start<Pos[i] & FACET.df$loc.end>Pos[i]] > (median.logr + std.logr), "loh", ".")
    }
    
}

muts <- cbind(muts, tcn, lcn, loh)

####################################################
################ Hotspot Annotation ################
####################################################
aa <- substr(muts$ANN....HGVS_P, 3, nchar(muts$ANN....HGVS_P))
aa.convert <- read.csv("/opt/src/runAbsolute/AminoAcidTable.csv", header = TRUE, stringsAsFactors = FALSE)
pos = first.aa.updated = last.aa.updated = '.'

for (k in 1:nrow(muts)){
    print(k)
    if(nchar(aa[k])>0){
        pos[k] = as.numeric(gsub("\\D", "", aa[k])) 
        first.aa = substr(aa, 1, 3) 
        first.aa.updated[k] = aa.convert$one.letter[first.aa[k] == aa.convert$three.letter]
        last.aa = substr(aa, nchar(aa)-2, nchar(aa))
        if (substr(last.aa[k], 3, 3) == "*"){last.aa.updated[k] = "*"}
        if (substr(last.aa[k], 2, 3) == "fs"){last.aa.updated[k] = "fs"}
        if (last.aa[k] %in% aa.convert$three.letter){last.aa.updated[k] = aa.convert$one.letter[last.aa[k] == aa.convert$three.letter]}
    }else{
        first.aa.updated[k] = "."
        pos [k] = ""
        last.aa.updated[k] = ""
    }
}

muts$AA <- paste(first.aa.updated, pos, last.aa.updated, sep = "")

## Annotate Hotspots
hotspots <- read.csv("/opt/src/runAbsolute/Hotspot_List_V3.csv", header = TRUE, stringsAsFactors = FALSE)

hotspot_id <- paste(hotspots$Gene, hotspots$AA, sep = "_")
muts$hs_ID = paste(muts$ANN....GENE, muts$AA, sep = "_")
muts$hotspot = ifelse(muts$hs_ID %in% hotspot_id, "TRUE", ".")
muts$ANN....EFFECT[muts$hotspot == "TRUE"] = paste(muts$ANN....EFFECT[muts$hotspot == "TRUE"], "hotspot", sep = "_")

write.csv(muts, "LM_results.csv", row.names = FALSE)

