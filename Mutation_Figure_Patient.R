## This will take the results file that was created using the combineMAF code and create 2 mutation files

## libraries
library(RColorBrewer)

## setwd('/Users/selenicp/Documents/Projects/Radiogenomics/WithinPatientSufam/ABS/')
results <- read.csv("LM_Results.csv", header = T, stringsAsFactors = F)

## If you are doing a sample set where there are multiple components you can use this to set the overall samples to get individual plots
## Otherwise this is just set to make all samples from one patient to make one plot.
results$Patient <- "A" #substr(results$NORMAL_SAMPLE, 1, nchar(results$NORMAL_SAMPLE)-1)
patients <- unique(results$Patient)


## If not all of the samples have mutations then this won't add those samples into the figure.
## Alternative method of getting sample list: 
## files <- list.files(pattern = "\\cncf.txt$")
## samples <- c(as.data.frame(strsplit(files, "_"))[1,])
## samples <- as.vector(unname(samples))
## samples <- as.character(unlist(samples))

for (b in 1:length(patients)){
    muts <- results[results$Patient == patients[b],]

    muts <- muts[!is.na(muts$Cancer_Cell_Fraction),]
    plot_file <- paste('Mutation_Heatmap', patients[b], ".pdf", sep = "")
    plot_file_CCF <- paste('CCF_Heatmap', patients[b], ".pdf", sep = "")

    
    muts$CCF_group <- "0"

    for (i in 1:nrow(muts)){
        if(muts$Cancer_Cell_Fraction[i] <= 1 & muts$Cancer_Cell_Fraction[i] > 0.8){muts$CCF_group[i] = 5}
        if(muts$Cancer_Cell_Fraction[i] <= 0.8 & muts$Cancer_Cell_Fraction[i] > 0.6){muts$CCF_group[i] = 4}
        if(muts$Cancer_Cell_Fraction[i] <= 0.6 & muts$Cancer_Cell_Fraction[i] > 0.4){muts$CCF_group[i] = 3}
        if(muts$Cancer_Cell_Fraction[i] <= 0.4 & muts$Cancer_Cell_Fraction[i] > 0.2){muts$CCF_group[i] = 2}
        if(muts$Cancer_Cell_Fraction[i] <= 0.2 & muts$Cancer_Cell_Fraction[i] > 0.05){muts$CCF_group[i] = 1}
    }
    
    
    muts <- muts[order(muts$ANN....GENE, decreasing = FALSE),]
    muts <- muts[!is.na(muts$Cancer_Cell_Fraction),]
    muts <- muts[order(muts$CCF_group, decreasing = TRUE),]
    muts <- muts[order(muts$TUMOR_SAMPLE, decreasing = FALSE),]
    
    ## Do any necessary filtering of mutation dataframe
                                        # muts <- subset(muts, subset = muts$patient == "RG6T") #| muts$TUMOR_SAMPLE == "AM2" | muts$TUMOR_SAMPLE == "AM3" | muts$TUMOR_SAMPLE == "AM4" | muts$TUMOR_SAMPLE == "AM5 - Primary_AME" | muts$TUMOR_SAMPLE == "AM6" | muts$TUMOR_SAMPLE == "AM7" | muts$TUMOR_SAMPLE == "AM32AME1")# & muts$TUMOR_SAMPLE != "AM8-Ips-Breast-Rec-Carc" & muts$TUMOR_SAMPLE != "AM5 - Axillary_Lymph_Node_Metastasis" & muts$TUMOR_SAMPLE != "AM46T2" & muts$TUMOR_SAMPLE != "AM5 - Primary_Carcarcinoma" & muts$TUMOR_SAMPLE != "AM8-LNE" & muts$TUMOR_SAMPLE != "AM8-LNM")
    muts <- subset(muts, subset = muts$ANN....EFFECT != "synonymous_variant") # & muts$ANN....EFFECT != "splice_region_variant&synonymous_variant" & muts$ANN....EFFECT != "downstream_ANN....GENE_variant" 
    ## & muts$ANN....EFFECT != "intron_variant" & muts$ANN....EFFECT != "frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant"
    ## & muts$ANN....EFFECT != "non_coding_exon_variant|synonymous_variant" & muts$ANN....EFFECT != "SYNONYMOUS_CODING" & muts$ANN....EFFECT != "splice_acceptor_variant&splice_region_variant&intron_variant" & muts$ANN....EFFECT != "Silent")
    ## muts = subset(muts, subset = muts$ANN....EFFECT != "downstream_gene_variant|synonymous_variant" & muts$ANN....EFFECT != "upstream_gene_variant|5_prime_UTR_variant")
    ## muts = subset(muts, subset = muts$lawrence == "TRUE" | muts$kandoth == "TRUE" | muts$cancer_gene_census == "TRUE")
    
    sample_names = as.list(sort(unique(muts$TUMOR_SAMPLE))) 
    if(length(sample_names) == 1){ sample_names = c(sample_names, "Spacer")}
    mutation_genes = unique(muts$ANN....GENE)
    if(length(mutation_genes) == 1){ mutation_genes = c(mutation_genes, "Spacer")}
    rownames(muts) = 1:nrow(muts)
    
    ## impact410_list = read.csv('/Users/selenicp/Documents/Code/IMPACT410_genes.csv', header = T, stringsAsFactors = F)
    ## impact341 = read.csv('/Users/burkek/Documents/Code/IMPACT341_genes.csv', header = TRUE, stringsAsFactors = F)
    ## impact410_list = unique(c(impact410_list$Approved.Symbol, impact410_list$HGNC.symbol))
    ## mutation_genes = unique(c(intersect(mutation_genes,impact410_list)))
    ## muts = muts[muts$ANN....GENE %in% mutation_genes,]
    
    return_gene_order = F
    sort_samples = T
    sort_genes = T
    show_sample_names = T
    TCGA = F
    remove_genes_with_no_mutation = F
    width = NULL 
    height = NULL
    include_percentages = T
    sample_name_col = "TUMOR_SAMPLE"
    
    ## ### Make a matrix of "blank" values the with nrow = #mutations and ncol = #samples
    mutation_heatmap <- matrix(9, nrow = sum(unlist(lapply(sample_names, length))), ncol = sum(unlist(lapply(mutation_genes, length))))
    rownames(mutation_heatmap) <- unlist(sample_names)
    colnames(mutation_heatmap) <- paste(mutation_genes)
    
    ## ### Make sure the sample and mutations are both in the list of gene mutations and gene samples
    if (!TCGA) { smallmaf <- muts[which(muts$ANN....GENE %in% unlist(mutation_genes) & muts$TUMOR_SAMPLE %in% unlist(sample_names)),]
    } else { 
        muts$id <- unlist(lapply(muts$TUMOR_SAMPLE, function(x){substr(x, 1, 12)}))
        print(head(muts$id))
        print(head(unlist(sample_names)))
        smallmaf <- muts[which(muts$Hugo_Symbol %in% unlist(mutation_genes) & muts$id %in% unlist(sample_names)),] 
    }
    
    ## ## Define categories for different Effects
    cat1 <- c("missense_variant_hotspot", "stop_gained_hotspot","missense_variant&splice_region_variant_hotspot")
    cat2 <- c("STOP_GAINED", "Nonsense_Mutation", "stop_gained&splice_region_variant", "stop_gained", "stop_gained&inframe_insertion", "stop_gained&disruptive_inframe_insertion", "stop_gained&disruptive_inframe_insertion&splice_region_variant", "stop_gained&inframe_insertion|stop_gained&inframe_insertion", "stop_gained&inframe_insertion&splice_region_variant", "stop_gained|stop_gained" )
    cat3 <- c("frameshift_indel","FRAME_SHIFT", "FRAME_SHIFT", "frameshift_variant&stop_lost", "frameshift_variant&stop_lost|frameshift_variant&stop_lost", "frameshift_variant&start_lost","frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant", "Frame_Shift_Del", "Frame_Shift_Ins", "frameshift_variant", "frameshift_variant|frameshift_variant", "frameshift_variant&stop_gained", "frameshift_variant&splice_region_variant", "frameshift_variant&splice_acceptor_variant&splice_region_variant&splice_region_variant&intron_variant", "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant", "frameshift_variant&stop_gained&splice_region_variant", "frameshift_variant&stop_gained|frameshift_variant&stop_gained","frameshift_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&splice_region_variant&splice_region_variant&intron_variant" )
    cat4 <- c("missense_snv_recurrent","missense_snv","NON_SYNONYMOUS_CODING", "Missense_Mutation", "missense_variant", "missense_variant&splice_region_variant", "missense_variant|missense_variant", "missense_variant|missense_variant|missense_variant|missense_variant|missense_variant|missense_variant|missense_variant|missense_variant|missense_variant", "missense_variant|missense_variant|missense_variant", "missense_variant&splice_region_variant|missense_variant&splice_region_variant")
    cat5 <- c("inframe_indel_recurrent","inframe_indel","CODON_CHANGE_PLUS_CODON_DELETION", "CODON_DELETION", "CODON_INSERTION", "In_Frame_Ins", "In_Frame_Del", "inframe_deletion|inframe_deletion", "disruptive_inframe_deletion", "disruptive_inframe_insertion", "inframe_deletion", "inframe_insertion", "disruptive_inframe_deletion&splice_region_variant", "inframe_deletion&splice_region_variant", "inframe_insertion&splice_region_variant", "disruptive_inframe_insertion|disruptive_inframe_insertion", "disruptive_inframe_insertion&splice_region_variant", "inframe_insertion|inframe_insertion", "disruptive_inframe_deletion|disruptive_inframe_deletion","conservative_inframe_deletion")
##     cat6 <- c("splice_donor_variant&intron_variant|splice_donor_variant&intron_variant", "splice_acceptor_variant&intron_variant|splice_acceptor_variant&intron_variant", "splice_acceptor_variant&splice_donor_variant&intron_variant", "splice_mut", "SPLICE_SITE_DONOR", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_REGION", "Splice_Site", "splice_donor_variant&intron_variant", "splice_acceptor_variant&intron_variant", "splicing", "splice_donor_variant&splice_region_variant&intron_variant", "splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&splice_region_variant&intron_variant", "splice_donor_variant&inframe_deletion&splice_region_variant&splice_region_variant&intron_variant", "splice_acceptor_variant&inframe_deletion&splice_region_variant&splice_region_variant&intron_variant", "Splice_Region", "splice_acceptor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|start_lost&inframe_deletion&splice_region_variant", "splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant","splice_acceptor_variant&splice_region_variant&intron_variant&non_coding_exon_variant", "splice_acceptor_variant&splice_region_variant&intron_variant")
    cat6 <- c("splice_region_variant&synonymous_variant|splice_region_variant&synonymous_variant","splice_donor_variant&intron_variant|splice_donor_variant&intron_variant","splice_acceptor_variant&intron_variant|splice_acceptor_variant&intron_variant","splice_acceptor_variant&splice_donor_variant&intron_variant","splice_mut","SPLICE_SITE_DONOR", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_REGION", "Splice_Site", "splice_donor_variant&intron_variant", "splice_acceptor_variant&intron_variant", "splicing", "splice_donor_variant&splice_region_variant&intron_variant",  "Splice_Region", "splice_acceptor_variant&splice_region_variant&intron_variant&non_coding_exon_variant", "splice_acceptor_variant&splice_region_variant&intron_variant","splice_region_variant&intron_variant")
    cat7 <- c("STOP_LOST", "START_LOST", "START_GAINED", "UTR_5_PRIME", "start_lost", "stop_lost", "start_lost&inframe_deletion", "stop_lost&splice_region_variant", "stop_lost&disruptive_inframe_insertion")
    cat8 <- c("upstream_gene_variant","downstream_gene_variant|synonymous_variant", "upstream_gene_variant|5_prime_UTR_variant","upstream_gene_variant|synonymous_variant")
    cat9 <- c("nonsense_snv","synonymous_variant","Silent", "splice_region_variant&synonymous_variant", "downstream_gene_variant", "intron_variant", "frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant", "non_coding_exon_variant|synonymous_variant", "SYNONYMOUS_CODING", "synonymous_variant|synonymous_variant", "splice_region_variant&synonymous_variant|splice_region_variant&non_coding_exon_variant", "intergenic_region", "intron_variant","intron_variant|downstream_gene_variant","intron_variant|intron_variant","intergenic_region|downstream_gene_variant","intron_variant|upstream_gene_variant")
    
    
    ## ### For each row read the Effect and create the type based on which category it fits in
    ## ### If there is an error because the mutation type is unknow just add it to the correct category above
    for (i in 1:nrow(smallmaf)) {
        if(!TCGA) { type = smallmaf$ANN....EFFECT[i] } else { type = smallmaf$Variant_Classification[i] }
        if (type %in% cat1) { type = 1
        } else if (type %in% cat2) { type = 2
        } else if (type %in% cat3) { type = 3
        } else if (type %in% cat4) { type = 4
        } else if (type %in% cat5) { type = 5 
        } else if (type %in% cat6) { type = 6
        } else if (type %in% cat7) { type = 7
        } else if (type %in% cat8) { type = 8
        } else if (type %in% cat9) { type = 9
        } else {print(paste(i,type,sep = "_"))
            stop("Mutation type not found")}
        print(paste(i,type,sep = "_"))
        
        if (!TCGA) { 
            if (mutation_heatmap[which(rownames(mutation_heatmap) == smallmaf$TUMOR_SAMPLE[i],), which(colnames(mutation_heatmap) == smallmaf$ANN....GENE[i])] > type) {
                mutation_heatmap[which(rownames(mutation_heatmap) == smallmaf$TUMOR_SAMPLE[i],), which(colnames(mutation_heatmap) == smallmaf$ANN....GENE[i])] = type
            } else { mutation_heatmap[which(rownames(mutation_heatmap) == smallmaf$TUMOR_SAMPLE[i]), which(colnames(mutation_heatmap) == smallmaf$Hugo[i])] = type }
        }
    }
    
    ## ##############################
    ## Order CCF-Group
    CCFgroup_heatmap <-matrix(0, nrow = sum(unlist(lapply(sample_names, length))), ncol = sum(unlist(lapply(mutation_genes, length))))
    rownames(CCFgroup_heatmap) <- unlist(sample_names)
    colnames(CCFgroup_heatmap) <- mutation_genes
    
    for (i in 1:nrow(smallmaf)) {
        type = smallmaf$CCF_group[i] 
        if (CCFgroup_heatmap[which(rownames(CCFgroup_heatmap) == smallmaf$TUMOR_SAMPLE[i],), which(colnames(CCFgroup_heatmap) == smallmaf$ANN....GENE[i])] < type) {
            CCFgroup_heatmap[which(rownames(CCFgroup_heatmap) == smallmaf$TUMOR_SAMPLE[i],), which(colnames(CCFgroup_heatmap) == smallmaf$ANN....GENE[i])] <- type}
    }

    i = 1
    for (i in 1 :nrow(mutation_heatmap)){
        mutation_heatmap = mutation_heatmap[, order(CCFgroup_heatmap[nrow(CCFgroup_heatmap)-(i-1),], decreasing = T) ]
        CCFgroup_heatmap = CCFgroup_heatmap[, order(CCFgroup_heatmap[nrow(CCFgroup_heatmap)-(i-1),], decreasing = T) ]
    }
    
    ## ######################################
    
    
    if (sort_samples) {
        mutation_heatmap <- do.call(rbind, lapply(sample_names, function(x) { m <- mutation_heatmap[which(rownames(mutation_heatmap) %in% x),, drop = FALSE]; m[do.call(order, transform(m)),]}))
    }
    rownames(mutation_heatmap) <- unlist(sample_names)
    
    ## ### Sort genes by the number of mutations that do not equal the blank category
    if(sort_genes) {
        print("Sorting genes")
        oo <- unlist(apply(mutation_heatmap,2, function(x){length(which(x != 9))}))
        print(oo)
        oo <- oo[which(!duplicated(names(oo)))]
        
        oo <- names(sort(oo, decr = TRUE))
        mutation_heatmap <- mutation_heatmap[,match(unlist(oo), colnames(mutation_heatmap))]	
        
    }
    
    
    if(remove_genes_with_no_mutation) { mutation_heatmap <- mutation_heatmap[,which(unlist(apply(mutation_heatmap,2, function(x){length(which(x != 10))})) != 0)] }
    order_heatmap = ifelse(mutation_heatmap == 9, 2, 1)
    i = 1
    for (i in 1 :ncol(mutation_heatmap)){
        mutation_heatmap = mutation_heatmap[ order(order_heatmap[,ncol(mutation_heatmap)-(i-1)]) , ]
        order_heatmap = order_heatmap[ order(order_heatmap[,ncol(mutation_heatmap)-(i-1)]) , ]
    }
    
    ## # 
    mutation_heatmap = t(mutation_heatmap)
    
    ## ### Choose color palette for table - set to colkey
    colkey <- cbind(1:9, c(brewer.pal(8, "Set1"), "white"))
    colkey[1,2] = "#CD1719" 
    colkey[2,2] = "#984EA3"
    colkey[3,2] = "#377EB8"
    colkey[4,2] = "#4DAF4A"
    colkey[5,2] = "#FF7F00"
    colkey[6,2] = "#FFFF33"
    colkey[7,2] = "#A65628"
    colkey[8,2] = "#808080"
    colkey[9,2] = "gray90"
    

    width = height = NULL
    if (is.null(width)) { width = 1+(length(unlist(mutation_genes))/1.5) }
    if (width<4){width = 4}
    if (is.null(height)) { height = 1+(length(unlist(sample_names)))/2 }
    if (height<4){height = 4}

    ## ### Create empty pdf
    pdf(plot_file, width = width, height = height)
    if (show_sample_names) { top = 8 } else {top = 2}
    if (include_percentages) { right = 4 } else { right = 1 }

    ## ### Plot figure
    par(oma = c(2,8,1,1), mar = c(2,5,top,right))
    image(mutation_heatmap, xaxt = "n", yaxt = "n", col = colkey[,2], zlim = c(1,9), xlab = "", ylab = "")
    axis(2, at = seq(0, 1, 1/(ncol(mutation_heatmap)-1)), labels = colnames(mutation_heatmap), las = 2, tick = FALSE, cex.axis = 2, font = 2, family = "sans")
    axis(3, at = seq(0, 1, 1/(nrow(mutation_heatmap)-1)), labels = rownames(mutation_heatmap), las = 2, tick = FALSE, cex.axis = 1.5, font = 3, family = "sans")
    abline(v = (0:nrow(mutation_heatmap)/(nrow(mutation_heatmap)-1)+(1/(2*(nrow(mutation_heatmap)-1)))), col = "white")
    abline(h = (0:ncol(mutation_heatmap)/(ncol(mutation_heatmap)-1)+(1/(2*(ncol(mutation_heatmap)-1)))), col = "white")
    dev.off()


    gene_order = rownames(mutation_heatmap)#[length(rownames(mutation_heatmap)):1]
    samp_order = colnames(mutation_heatmap)

    if(return_gene_order) { list(rev(colnames(mutation_heatmap))) }



    return_gene_order = F
    sort_samples = T
    sort_genes = T
    show_sample_names = T
    TCGA = F
    remove_genes_with_no_mutation = F
    ## width = NULL 
    ## height = NULL
    include_percentages = F
    sample_name_col = "TUMOR_SAMPLE"
    mutation_genes = unique(muts$ANN....GENE)

    ## ### Make a matrix of "blank" values the with nrow = #mutations and ncol = #samples
    mutation_heatmap <- matrix(1, nrow = sum(unlist(lapply(samp_order, length))), ncol = sum(unlist(lapply(gene_order, length))))
    rownames(mutation_heatmap) <- samp_order
    colnames(mutation_heatmap) <- gene_order

    ## ### Make sure the sample and mutations are both in the list of gene mutations and gene samples
    if (!TCGA) { smallmaf <- muts[which(muts$ANN....GENE %in% unlist(mutation_genes) & muts$TUMOR_SAMPLE %in% unlist(sample_names)),]
    }else { 
        muts$id <- unlist(lapply(muts$TUMOR_SAMPLE, function(x){substr(x, 1, 12)}))
        print(head(muts$id))
        print(head(unlist(sample_names)))
        smallmaf <- muts[which(muts$Hugo_Symbol %in% unlist(mutation_genes) & muts$id %in% unlist(sample_names)),] 
    }

    ## ### For each row read the Effect and create the type based on which category it fits in
    for (i in 1:nrow(smallmaf)) {
        if(!TCGA) { type = smallmaf$Cancer_Cell_Fraction[i] } else { type = smallmaf$Variant_Classification[i] }
        if (type == 0) { type = 1
        } else if (type >0 & type<= 0.05) { type = 2
        } else if (type >0.05 & type<= 0.2) { type = 3
        } else if (type >0.2 & type<= 0.4) { type = 4
        } else if (type >0.4 & type<= 0.6) { type = 5 
        } else if (type >0.6 & type<= 0.8) { type = 6
        } else if (type >0.8 & type<= 1) { type = 7
        } else { print("Mutation type not found") }
        print(paste(i,type,sep = "_"))
        
        if (mutation_heatmap[which(rownames(mutation_heatmap) == smallmaf$TUMOR_SAMPLE[i],), which(colnames(mutation_heatmap) == smallmaf$ANN....GENE[i])] < type) {
            mutation_heatmap[which(rownames(mutation_heatmap) == smallmaf$TUMOR_SAMPLE[i],), which(colnames(mutation_heatmap) == smallmaf$ANN....GENE[i])] <- type}
    }

    if(remove_genes_with_no_mutation) { mutation_heatmap <- mutation_heatmap[,which(unlist(apply(mutation_heatmap,2, function(x){length(which(x != 10))})) != 0)] }

    df = data.frame(matrix(ncol = 0, nrow = nrow(mutation_heatmap)))
    for (i in 1:length(gene_order)){
        print(i)
        df = cbind.data.frame(df, mutation_heatmap[ , colnames(mutation_heatmap) == gene_order[i]])
        colnames(df)[i] = gene_order[i] 
    }
    rownames(df) = rownames(mutation_heatmap)
    newdf = data.frame(matrix(ncol = 0, nrow = ncol(mutation_heatmap)))
    for (i in 1:length(samp_order)){
        print(i)
        newdf = rbind.data.frame(newdf, df[rownames(mutation_heatmap) == samp_order[i],])
        rownames(newdf)[i] = samp_order[i] 
    }
    colnames(newdf) = colnames(df)
    mutation_heatmap <- as.matrix(newdf)


    ## Choose color palette for table - set to colkey
    colkey <- cbind(1:7, c("grey90", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5", "#08519C", "#08306B"))

    mutation_heatmap = t(mutation_heatmap)


    ## Make a matrix of "blank" values the with nrow = #mutations and ncol = #samples
    LOH_heatmap <- matrix(0, ncol = sum(unlist(lapply(gene_order, length))), nrow = sum(unlist(lapply(samp_order, length))))
    colnames(LOH_heatmap) <- gene_order
    rownames(LOH_heatmap) <- samp_order

    for (i in 1:nrow(smallmaf)) {
        type = smallmaf$loh[i]
        if (type == "loh") { type = 1
                                        # } else if (type >0 & type<= 0.05) { type = 2
        } else { type = 0 }
        print(paste(i,type,sep = "_"))
        
        if (LOH_heatmap[which(rownames(LOH_heatmap) == smallmaf$TUMOR_SAMPLE[i]), which(colnames(LOH_heatmap) == smallmaf$ANN....GENE[i],)] < type) {
            LOH_heatmap[which(rownames(LOH_heatmap) == smallmaf$TUMOR_SAMPLE[i]), which(colnames(LOH_heatmap) == smallmaf$ANN....GENE[i],)] <- type}
    }

    ## Make a matrix of "blank" values the with nrow = #mutations and ncol = #samples
    Clonal_heatmap <- matrix(0, ncol = sum(unlist(lapply(gene_order, length))), nrow = sum(unlist(lapply(samp_order, length))))
    colnames(Clonal_heatmap) <- gene_order 
    rownames(Clonal_heatmap) <- samp_order

    for (i in 1:nrow(smallmaf)) {
                                        # type = smallmaf$clonality[i]
                                        # if (type == "clonal") { type = 1
        type = smallmaf$Clonal_Status[i]
        if (type == "Clonal") { type = 1
                                        # } else if (type >0 & type <= 0.05) { type = 2
        } else { type = 0 }
        print(paste(i,type,sep = "_"))
        
        if (Clonal_heatmap[which(rownames(Clonal_heatmap) == smallmaf$TUMOR_SAMPLE[i]), which(colnames(Clonal_heatmap) == smallmaf$ANN....GENE[i],)] < type) {
            Clonal_heatmap[which(rownames(Clonal_heatmap) == smallmaf$TUMOR_SAMPLE[i]), which(colnames(Clonal_heatmap) == smallmaf$ANN....GENE[i],)] <- type}
    }

    ## Create empty pdf
    pdf(plot_file_CCF, height = height, width = width)
    if (show_sample_names) { top = 8 } else {top = 2}
    if (include_percentages) { right = 4 } else { right = 1 }



    ## Clonal_heatmap <- Clonal_heatmap[,rowSums(mutation_heatmap != 1) > 1]
    ## LOH_heatmap <- LOH_heatmap[,rowSums(mutation_heatmap != 1) > 1]
    ## mutation_heatmap <- mutation_heatmap[rowSums(mutation_heatmap != 1) > 1,]

    ## Create empty pdf
    ## width = height = NULL
    ## if (is.null(width)) { width = 1+nrow(mutation_heatmap)/1.5 }
    ## if (width<4){width = 4}
    ## if (is.null(height)) { height = 1+(length(unlist(sample_names)))/2 }
    ## if (height<4){height = 4}

    ## pdf(plot_file_CCF, height = height, width = width)
    ## if (show_sample_names) { top = 8 } else {top = 2}
    ## if (include_percentages) { right = 4 } else { right = 1 }



    ## ### Plot figure
    par(oma = c(2,8,1,1), mar = c(2,5,top,right))
    ## par(mar = c(2,4,2,2))
    image(mutation_heatmap, xaxt = "n", yaxt = "n", col = colkey[,2], zlim = c(1,7), xlab = "", ylab = "")
    ## axis(1, at = seq(0, 1, 1/(ncol(mutation_heatmap)-1)), labels = colnames(mutation_heatmap), las = 2, tick = FALSE)
    axis(2, at = seq(0, 1, 1/(ncol(mutation_heatmap)-1)), labels = colnames(mutation_heatmap), las = 2, tick = FALSE, cex.axis = 2, font = 3, family = "sans")
    axis(3, at = seq(0, 1, 1/(nrow(mutation_heatmap)-1)), labels = rownames(mutation_heatmap), las = 2, tick = FALSE, cex.axis = 1.2, font = 2, family = "sans")
    if(include_percentages) {
        axis(4, at = seq(0, 1, 1/(ncol(mutation_heatmap)-1)), labels = unlist(apply(mutation_heatmap, 2, function(x) { paste(round(100*length(which(x != 1))/length(x), 0), "%")})), las = 2, tick = FALSE)}
    ## abline(v = (0:nrow(mutation_heatmap)/(nrow(mutation_heatmap)-1)+(1/(2*(nrow(mutation_heatmap)-1)))), col = "white")
    ## abline(h = (0:ncol(mutation_heatmap)/(ncol(mutation_heatmap)-1)+(1/(2*(ncol(mutation_heatmap)-1)))), col = "white")

    row_pos = 0:nrow(mutation_heatmap)/(nrow(mutation_heatmap)-1)
    row_pos = row_pos[!row_pos>1]
    col_pos = 0:ncol(mutation_heatmap)/(ncol(mutation_heatmap)-1)
    col_pos = col_pos[!col_pos>1]
    ## LOH_heatmap = LOH_heatmap[,ncol(LOH_heatmap):1]
    ## LOH_heatmap = t(LOH_heatmap)
    LOH_heatmap = LOH_heatmap[nrow(LOH_heatmap):1,]
    for (i in 1:nrow(LOH_heatmap)){
        row_plot = row_pos[as.vector(LOH_heatmap[i,] == 1)]
        points( row_plot, rep(col_pos[(length(samp_order)+1)-i],length(row_plot)), pch = "/", col = "white", cex = 3)
    }
    Clonal_heatmap = Clonal_heatmap[nrow(Clonal_heatmap):1,]
    for (i in 1:nrow(Clonal_heatmap)){
        row_plot = row_pos[as.vector(Clonal_heatmap[i,] == 1)]
        points( row_plot, rep(col_pos[(length(samp_order)+1)-i],length(row_plot)), pch = 16, col = "yellow", cex = 2)
    }

    dev.off()
}
