#Process the recent and recurrent selection summary statistics -- FINAL
#Using DHEW and only considering the results from DH program (TajD, nFWH, EW, HEW, and DHEW)
#September 5, 2017 and edited extensively on Nov 27th
#Joo Hyun Im (ji72)

########################################################################################################
#A. Get the data
########################################################################################################

#Delete previous data
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/")

#Indicate which species we're working with:
#indicator = "Dmel"
indicator = "Dsim"


#Get the DH data
get_the_data = function(indicator){
    
    #Read-in the corresponding data
    if (indicator == "Dmel"){
        result.DHEW = read.table("final_data/Dmel_DHEW_extracted_compiled_Nov_2017_cleaned.txt", sep="\t", header=T, stringsAsFactors = F) # genes (the genes without enough control genes have not been excluded)
        
    }
    else{
        result.DHEW = read.table("final_data/Dsim_DHEW_extracted_compiled_Nov_2017_cleaned.txt", sep="\t", header=T, stringsAsFactors = F) #477 genes (the genes without enough control genes have not been excluded)
    }
    
    #Add functional categories
    if (indicator == "Dmel"){
        result.label = read.table("final_data/internalization_list_Dmel_Nov_2017_final.txt", sep="\t", header=T, stringsAsFactors = F) # genes (the genes without enough control genes have not been excluded)
        
    }
    else{
        result.label = read.table("final_data/internalization_list_Dsim_Nov_2017_final.txt", sep="\t", header=T, stringsAsFactors = F) #477 genes (the genes without enough control genes have not been excluded)
    }
    
    result.total =c()
    for (i in 1:dim(result.DHEW)[1]){
        gene_id = as.character(result.DHEW[i,1])
        if (gene_id %in% result.label[,1]){ #if the gene is not filtered out due to a lack/smaller # of control genes
            result.thread = c(as.character(result.label[which(result.label$gene_id == gene_id),c(1:2)]), as.numeric(result.DHEW[i,2:12]), as.character(result.label[which(result.label$gene_id == gene_id),c(3:6)]))
            result.total = rbind(result.total, result.thread)
        }
        else{
            #cat("missing genes: ", gene_id, "\n")
        }
    }
    
    colnames(result.total) = c("gene_id", "gene_name", "S", "thetaW", "TajD", "TajD_pval", "nFWH", "nFWH_pval", "EW", "EW_pval", "DH_pval", "HEW_pval", "DHEW_pval", "target.control.status","type","subtype","matching.target")
    result.total = data.frame(result.total)
    
    result.total$S = as.numeric(levels(result.total$S))[result.total$S]
    result.total$thetaW = as.numeric(levels(result.total$thetaW))[result.total$thetaW]
    result.total$TajD = as.numeric(levels(result.total$TajD))[result.total$TajD]
    result.total$TajD_pval = as.numeric(levels(result.total$TajD_pval))[result.total$TajD_pval]
    result.total$nFWH = as.numeric(levels(result.total$nFWH))[result.total$nFWH]
    result.total$nFWH_pval = as.numeric(levels(result.total$nFWH_pval))[result.total$nFWH_pval]
    result.total$EW = as.numeric(levels(result.total$EW))[result.total$EW]
    result.total$EW_pval = as.numeric(levels(result.total$EW_pval))[result.total$EW_pval]
    result.total$DH_pval = as.numeric(levels(result.total$DH_pval))[result.total$DH_pval]
    result.total$HEW_pval = as.numeric(levels(result.total$HEW_pval))[result.total$HEW_pval]
    result.total$DHEW_pval = as.numeric(levels(result.total$DHEW_pval))[result.total$DHEW_pval]
    
    result.total = result.total[order(result.total$target.control.status, decreasing = T),]
    return(result.total)
}
result.total = get_the_data(indicator)
#write.table(result.total, file="final_data/Dmel_DHEW_extracted_compiled_Nov_2017_cleaned_with_function.txt", quote=F, row.names = F, col.names = T, sep="\t")
write.table(result.total, file="final_data/Dsim_DHEW_extracted_compiled_Nov_2017_cleaned_with_function.txt", quote=F, row.names = F, col.names = T, sep="\t")


#Check the number of genes in the data
cat("Count the number of genes in the following species: ", indicator, "\n")
result.target = result.total[which(result.total$target.control.status == "target"),]
result.cont = result.total[which(result.total$target.control.status == "control"),]
cat("The number of target genes for ", indicator, " is: ", dim(result.target)[1], "\n")
cat("The number of control genes for ", indicator, " is: ", dim(result.cont)[1], "\n")
cat("Check if each target gene has at least 3 control genes")
table(result.cont$matching.target)

#Check if there are control genes for each target gene
result.cont$matching.target = factor(result.cont$matching.target); matching.target.genes = levels(result.cont$matching.target)
dim(result.target)[1]; length(matching.target.genes) #these should be the same


#Get the MK data
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/final_data/")

###I. Basic processing
#indicator = "Dmel"
indicator = "Dsim"

basic_processing = function(indicator){
    
    #1. Before reading in the file, I changed the file names listed in the document below to just have the FBgn names.
    #2. Read in the data
    read_data <- function(indicator){
        if (indicator == "Dmel"){
            return(read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/final_data/Dmels_Nov_2017_MK_results_cleaned.txt", header=F, stringsAsFactors = F))
        } else if (indicator == "Dsim"){
            return(read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/final_data/Dsims_Nov_2017_MK_results_cleaned.txt", header=F, stringsAsFactors = F))
        }
    }
    mk_test_results = read_data(indicator)
    colnames(mk_test_results) = c("gene_id","Pn", "Dn", "Ps", "Ds","MKcodons","FETpval")
    number_of_genes = dim(mk_test_results)[1]
    
    
    #3. Add the gene and function names
    read_func_data <- function(indicator){
        if (indicator == "Dmel"){
            return(read.table("internalization_list_Dmel_Nov_2017_final.txt", header=T, sep ='\t', stringsAsFactors = F))
        } else if (indicator == "Dsim"){
            return(read.table("internalization_list_Dsim_Nov_2017_final.txt", header=T, sep='\t', stringsAsFactors = F))
        }
    }
    gene_list = read_func_data(indicator)
    
    result_total =c(); result_thread =c()
    add_gene_name_and_function = function(){
        for (i in 1:dim(mk_test_results)[1]){
            gene_id = as.character(mk_test_results[i,1])
            if (gene_id %in% gene_list[,1]){
                result_thread = c(as.character(gene_list[which(gene_list$gene_id == gene_id),c(1:2)]), as.numeric(mk_test_results[i,2:7]), as.character(gene_list[which(gene_list$gene_id == gene_id),c(3:6)]))
                result_total = rbind(result_total, result_thread)
            }
        }
        return(result_total)
    }
    mk_test_results_final = add_gene_name_and_function()
    colnames(mk_test_results_final) = c("gene_id","gene_name","Pn", "Dn", "Ps", "Ds", "MKcodons","FETpval", "target_control_status", "type", "subtype","matching_target")
    mk_test_results_final[,c(3:8)] = as.numeric(mk_test_results_final[,c(3:8)])
    mk_test_results_final = data.frame(mk_test_results_final)
    
    
    #4.Multiple testing correction
    mk_test_results_final$FETpval = as.numeric(as.character(mk_test_results_final$FETpval))
    corrected_pval = p.adjust(mk_test_results_final$FETpval, method = "BH")
    mk_test_results_final_with_corrected_pval = cbind(mk_test_results_final, corrected_pval)
    
    
    #5. Calculate DoS (DoS = Dn/(Dn + Ds) − Pn/(Pn + Ps))
    mk_test_results_final_with_corrected_pval$Dn = as.numeric(as.character(mk_test_results_final_with_corrected_pval$Dn))
    mk_test_results_final_with_corrected_pval$Pn = as.numeric(as.character(mk_test_results_final_with_corrected_pval$Pn))
    mk_test_results_final_with_corrected_pval$Ds = as.numeric(as.character(mk_test_results_final_with_corrected_pval$Ds))
    mk_test_results_final_with_corrected_pval$Ps = as.numeric(as.character(mk_test_results_final_with_corrected_pval$Ps))
    calculate_DoS = function(data){
        DoS_first = (data$Dn) / (data$Dn + data$Ds); DoS_second = (data$Pn) / (data$Pn + data$Ps)
        DoS = DoS_first - DoS_second
        data_with_DoS = cbind(data, DoS)
        
        return(data_with_DoS[order(data_with_DoS$corrected_pval), ])
    }
    mk_test_results_final_with_corrected_pval_and_DoS = calculate_DoS(mk_test_results_final_with_corrected_pval)
    
    
    #6. Select the ones with significant p-values
    ordered_and_selected = mk_test_results_final_with_corrected_pval_and_DoS[which(mk_test_results_final_with_corrected_pval_and_DoS$corrected_pval < 0.05),]
    #print(ordered_and_selected)
    #results: Rel was the only one with the significant p-value.
    
    #7. Save the data
    return(mk_test_results_final_with_corrected_pval_and_DoS)
    
}
result.total = basic_processing(indicator)

#Check the number of genes in the data
cat("Count the number of genes in the following species: ", indicator, "\n")
result.target = result.total[which(result.total$target_control_status == "target"),]
result.cont = result.total[which(result.total$target_control_status == "control"),]
cat("The number of target genes for ", indicator, " is: ", dim(result.target)[1], "\n")
cat("The number of control genes for ", indicator, " is: ", dim(result.cont)[1], "\n")

#Save the results
#write.table(result.total, file="Dmels_Nov_2017_MK_results_analyzed_final.txt", quote=F, col.names = T, row.names = F, sep="\t")
write.table(result.total, file="Dsims_Nov_2017_MK_results_analyzed_final.txt", quote=F, col.names = T, row.names = F, sep="\t")

#Check if there are control genes for each target gene
result.cont$matching.target = factor(result.cont$matching_target); matching.target.genes = levels(result.cont$matching_target)
dim(result.target)[1]; length(matching.target.genes) #these should be the same
setdiff(matching.target.genes, as.character(result.target$gene_id)) 

#Check if each target gene has at least 3 control genes
table(result.cont$matching_target) 

#Dsim has the following problematic genes. I've manually fixed them in Dsims_Nov_2017_MK_results_analyzed_final.txt:
#"FBgn0011274" "FBgn0020377" "FBgn0022131" "FBgn0028741" "FBgn0030926" "FBgn0043577", #"FBgn0027594", "FBgn0029943", "FBgn0030926", "FBgn0035975"



########################################################################################################
#B. Comparison of statistics between Internalization, Recognition, and Signaling
########################################################################################################

#Francoise Vermelyan's suggestion:
#Calculate median(control) - target for a given target gene and perform a 1-sample t-test to see if the MEAN of these numbers is significantly different from 0

#Input
result.total.r1 = result.total #thetwW, TajD, nFWH, EW, #save it as a separate dataset
#result.total.r1 = read.table(paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/final_data/",indicator,"s_Nov_2017_MK_results_analyzed_final.txt", sep=""), header=T, sep='\t') ##DoS

#Put together autophagy and phagocytosis genes
result.total.r1$type = revalue(result.total.r1$type, c("Autophagy" = "Internalization", "Phagocytosis" = "Internalization", "Both" = "Internalization")) #change [category] into 'Internalization'
result.total.r1$type = factor(result.total.r1$type, levels=unique(result.total.r1$type))

one_sample_t_test = function(bio_pathway, test_choice){
    cat("Perform a 1-sample t-test for the ", test_choice, " statistic in ", bio_pathway, " genes.", "\n")
    
    list.of.target.genes = levels(result.total.r1$matching.target) #list of case genes that have control genes (we got the names from matching case column)
    diff.table = data.frame()
    matching_target_number =17
    
    if(test_choice == "TajD"){
        test_number = 5 #TajD
    }
    else if (test_choice == "nFWH"){
        test_number = 7
    }
    else if (test_choice == "thetaW"){
        test_number = 4 
    }
    else if (test_choice == "EW"){
        test_number = 9 
    }
    else if (test_choice == "DoS"){
        test_number = 14
        list.of.target.genes = levels(result.total.r1$matching_target)
        matching_target_number =12
    }
    
    for (n in 1:length(list.of.target.genes)){
        pattern = list.of.target.genes[n]
        data.for.pattern.target = result.total.r1[grep(pattern, result.total.r1[,1]),]
        data.for.pattern.target.test = data.for.pattern.target[,test_number]
        data.for.pattern.target.name = data.for.pattern.target$gene_name
        data.for.pattern.cont = result.total.r1[grep(pattern, result.total.r1[,matching_target_number]),] #picking the matching target
        data.for.pattern.cont.test.median = median(data.for.pattern.cont[,test_number], na.rm = T)
        type = data.for.pattern.cont$type[1]; subtype = data.for.pattern.cont$subtype[1]
        diff.median = (data.for.pattern.target.test - data.for.pattern.cont.test.median) # difference between target
        diff.table[n,1] = pattern; diff.table[n,2] = data.for.pattern.target.name; diff.table[n,3] = diff.median; diff.table[n,4] = type; diff.table[n,5] = subtype
    } #case minus control
    colnames(diff.table) = c("target.id", "target.name", "diff.median", "type", "subtype") #89 target genes
    diff.table.sorted = diff.table[order(diff.table$diff.median, decreasing=T),]
    diff.table.function = diff.table.sorted[which(diff.table.sorted$type == bio_pathway),]
    cat("Mean: ", mean(diff.table.function$diff.median, na.rm = T), "\n")
    cat("Number of genes: ", dim(diff.table.function)[1], "\n")
    
    #Create a confidence interval
    #Check the normality assumption
    cat("P-value for checking the normality assumption (Shapiro): ", shapiro.test(c(diff.table.function$diff.median))$p.value, "\n")
    if (shapiro.test(c(diff.table.function$diff.median))$p.value < 0.05){
        cat("Normality: not normally distributed", "\n") #If <0.05, NOT NORMALLY DISTRIBUTED
    }
    else {
        cat("Normality: normally distributed", "\n")
    }
    
    t_test_result = t.test(diff.table.function$diff.median, mu=0)
    median_list = na.omit(diff.table.function$diff.median)
    hist(median_list, main = paste("Hist of ", test_choice, sep = ""))
    
    return(t_test_result) #The intervals calculated using the basic bootstrap method.

    
}

par(mfrow=c(1,4))

#Internalization
one_sample_t_test("Internalization", "thetaW")
one_sample_t_test("Internalization", "TajD")
one_sample_t_test("Internalization", "nFWH")
one_sample_t_test("Internalization", "EW")
one_sample_t_test("Internalization", "DoS")

#Signaling
one_sample_t_test("Signaling", "thetaW")
one_sample_t_test("Signaling", "TajD")
one_sample_t_test("Signaling", "nFWH")
one_sample_t_test("Signaling", "EW")
one_sample_t_test("Signaling", "DoS")

#Recognition
one_sample_t_test("Recognition", "thetaW")
one_sample_t_test("Recognition", "TajD")
one_sample_t_test("Recognition", "nFWH")
one_sample_t_test("Recognition", "EW")
one_sample_t_test("Recognition", "DoS")



########################################################################################################
#C. Create comprehensive data file for recent selection
########################################################################################################

rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/")

#Indicate which species we're working with:
#indicator = "Dmel"
indicator = "Dsim"

result.total = read.table(paste("final_data/", indicator, "_DHEW_extracted_compiled_Nov_2017_cleaned_with_function.txt", sep=""), header=T, sep='\t')
colnames(result.total)[14] = c("target_control_status") #unify the colname


#1) Calculate the adjusted p-values 

#only look at the internalization genes)
result.target = result.total[which(result.total$target_control_status == "target"),]
result.target.simplified = result.target[which(result.target$type %in% c("Autophagy", "Phagocytosis", "Both")),c(1,2,4,5,6,7,8,9,10,11,12,13,15,16)]

#Calculate the adjusted p-values
result.target.simplified[,15] = p.adjust(result.target.simplified$TajD_pval, method = "BH")
result.target.simplified[,16] = p.adjust(result.target.simplified$nFWH_pval, method = "BH")
result.target.simplified[,17] = p.adjust(result.target.simplified$EW_pval, method = "BH")
result.target.simplified[,18] = p.adjust(result.target.simplified$DH_pval, method = "BH")
result.target.simplified[,19] = p.adjust(result.target.simplified$HEW_pval, method = "BH")
result.target.simplified[,20] = p.adjust(result.target.simplified$DHEW_pval, method = "BH")
colnames(result.target.simplified)[15:20] = c("TajD_pval_adj", "nFWH_pval_adj", "EW_pval_adj", "DH_pval_adj", "HEW_pval_adj", "DHEW_pval_adj")

#Save the p-values
write.table(result.target.simplified, file = paste("final_data/", indicator, "_DHEW_extracted_compiled_Nov_2017_cleaned_with_function_adj_pval.txt", sep=""), quote=F, row.names = F, col.names = T)



#2) Calculate the difference between median of control genes and target gene and rank them (X out of 74).

result.target = result.total[which(result.total$target_control_status == "target"),]
result.cont = result.total[which(result.total$target_control_status == "control"),]
result.target.internalized = result.target[which(result.target$type %in% c("Autophagy", "Phagocytosis", "Both")),] #only get internalization genes
list.of.target.genes = droplevels(result.target.internalized$gene_id)

rank_the_difference = function(list.of.target.genes, statistic){
    diff.table = data.frame()
    for (n in 1:length(list.of.target.genes)){
        pattern = list.of.target.genes[n]
        data.for.pattern.target = result.total[grep(pattern, result.total[,1]),]
        
        if (statistic == "TajD"){
            data.for.pattern.target.statistic = data.for.pattern.target$TajD
            data.for.pattern.target.statistic.pval = data.for.pattern.target$TajD_pval
            
        }
        else if (statistic == "nFWH"){
            data.for.pattern.target.statistic = data.for.pattern.target$nFWH
            data.for.pattern.target.statistic.pval = data.for.pattern.target$nFWH_pval
        }
        else if (statistic == "EW"){
            data.for.pattern.target.statistic = data.for.pattern.target$EW
            data.for.pattern.target.statistic.pval = data.for.pattern.target$EW_pval
        }
        
        data.for.pattern.target.name = data.for.pattern.target$gene_name
        data.for.pattern.cont = result.total[grep(pattern, result.total[,17]),] #matching.target
        
        if (statistic == "TajD"){
            data.for.pattern.cont.statistic.median = median(data.for.pattern.cont$TajD)
        }
        else if (statistic == "nFWH"){
            data.for.pattern.cont.statistic.median = median(data.for.pattern.cont$nFWH)
        }
        else if (statistic == "EW"){
            data.for.pattern.cont.statistic.median = median(data.for.pattern.cont$EW)
        }
        
        type = data.for.pattern.cont$type[1]; subtype = data.for.pattern.cont$subtype[1]
        diff.median = (data.for.pattern.target.statistic - data.for.pattern.cont.statistic.median) # difference between target
        
        diff.table[n,1] = pattern; diff.table[n,2] = data.for.pattern.target.name
        diff.table[n,3] = data.for.pattern.target.statistic
        diff.table[n,4] = data.for.pattern.cont.statistic.median
        diff.table[n,5] = diff.median; diff.table[n,6] = type; diff.table[n,7] = subtype
        diff.table[n,8] = data.for.pattern.target.statistic.pval
        
    } #target minus control
    
    colnames(diff.table) = c("target_id", "target_name", "target_value", "med_cont_value", "diff_median", "type", "subtype", "target_pval")
    diff.table.sorted = diff.table[order(diff.table$diff_median, decreasing=F),]
    percentage_ranking = c(1:length(list.of.target.genes))/length(list.of.target.genes)*100
    diff.table.sorted = cbind(diff.table.sorted, c(1:length(list.of.target.genes)), percentage_ranking)
    colnames(diff.table.sorted)[9] = "ranking"; colnames(diff.table.sorted)[10] = "percentage_ranking"
    
    return(diff.table.sorted)
}
ranking_of_median_difference_TajD = rank_the_difference(list.of.target.genes, "TajD")
ranking_of_median_difference_nFWH = rank_the_difference(list.of.target.genes, "nFWH")
ranking_of_median_difference_EW = rank_the_difference(list.of.target.genes, "EW")



#3) Create a distribution of [gene 1 - median(rest of genes)] and see where my [target - median(control)] lands — breaks the pairing/grouping of focal vs control, but you can actually create a distribution (X out of 375)

result.target.internalized = result.target[which(result.target$type %in% c("Autophagy", "Phagocytosis", "Both")),] #only get internalization genes
list.of.target.genes = droplevels(result.target.internalized$gene_id)
result.total.internalized = result.total[which(result.total$type %in% c("Autophagy", "Phagocytosis", "Both")),] #only get internalization genes with their controls


#Create a null distribution
create_a_null_distribution = function(list.of.target.genes, statistic){
    data.for.dist = result.total.internalized #internalization genes only + their controls
    data.for.dist.target.only = data.for.dist[which(data.for.dist$target_control_status == "target"),]
    data.for.dist.cont.only = data.for.dist[which(data.for.dist$target_control_status == "control"),]
    data.for.dist.target.only$matching.target = data.for.dist.target.only$gene_id
    data.for.dist = rbind(data.for.dist.target.only, data.for.dist.cont.only)
    
    data.for.dist$matching.target = factor(data.for.dist$matching.target)
    list.of.target.genes = levels(data.for.dist$matching.target)
    
    if (statistic == "TajD"){
        column_number = 5
    }
    else if (statistic == "nFWH"){
        column_number = 7
    }
    else if (statistic == "EW"){
        column_number = 9
    }
    
    median.table = c()
    for (n in 1:length(list.of.target.genes)){
        cluster = data.for.dist[which(data.for.dist$matching.target == list.of.target.genes[n]),]
        number.of.genes.in.cluster = dim(cluster)[1]
        for (o in 2:number.of.genes.in.cluster){ #remove the true [target - median(control)] case
            target.data = cluster[o,column_number]; median.control.data = median(cluster[-o,column_number])
            diff.median = (target.data - median.control.data)
            median.table = rbind(median.table, diff.median)
        }
    }
    median.table.ordered = median.table[order(median.table[,1], decreasing=T),] #target minus control
    print(shapiro.test(as.numeric(unlist(median.table.ordered))))
    hist(median.table.ordered, main=paste("Null distribution of [gene 1 - median(rest of genes)]: ", statistic, sep=""))
    
    return(median.table.ordered)
}
par(mfrow=c(3,1))
null_distribution_TajD = create_a_null_distribution(list.of.target.genes, "TajD") #240
null_distribution_nFWH = create_a_null_distribution(list.of.target.genes, "nFWH") #240
null_distribution_EW = create_a_null_distribution(list.of.target.genes, "EW") #240

#Compare the true [target - median(control)] to the null distribution
compare_to_null = function(list.of.target.genes, statistic, null.dist){
    data.for.dist = result.total.internalized
    data.for.dist.target.only = data.for.dist[which(data.for.dist$target_control_status == "target"),]
    data.for.dist.cont.only = data.for.dist[which(data.for.dist$target_control_status == "control"),]
    data.for.dist.target.only$matching.target = data.for.dist.target.only$gene_id
    data.for.dist = rbind(data.for.dist.target.only, data.for.dist.cont.only)
    
    data.for.dist$matching.target = factor(data.for.dist$matching.target)
    list.of.target.genes = levels(data.for.dist$matching.target)
    
    if (statistic == "TajD"){
        column_number = 5
    }
    else if (statistic == "nFWH"){
        column_number = 7
    }
    else if (statistic == "EW"){
        column_number = 9
    }
    
    median.table = c()
    for (n in 1:length(list.of.target.genes)){
        cluster = data.for.dist[which(data.for.dist$matching.target == list.of.target.genes[n]),]
        number.of.genes.in.cluster = dim(cluster)[1]
        target.data = cluster[1,column_number]; median.control.data = median(cluster[-1,column_number])
        diff.median = (target.data - median.control.data)
        
        #find out where the true difference stands in the null distribution
        null.dist.with.true = c(null.dist, diff.median)
        ranking = rank(null.dist.with.true)[length(null.dist)+1] #ranking of the true difference - the last one
        ranking_percentage = c(ranking/length(null.dist))*100
        gene_name = as.character(data.for.dist[which(data.for.dist$gene_id == list.of.target.genes[n]),2])
        type = as.character(data.for.dist[which(data.for.dist$gene_id == list.of.target.genes[n]),15])
        subtype = as.character(data.for.dist[which(data.for.dist$gene_id == list.of.target.genes[n]),16])
        diff.median.and.rank = c(list.of.target.genes[n], gene_name, type, subtype, diff.median, ranking, ranking_percentage)
        
        median.table = rbind(median.table, diff.median.and.rank)
    }
    colnames(median.table) = c("gene_id", "gene_name", "type", "subtype", "True_median_difference", "Ranking", "percentage_ranking")
    median.table = median.table[order(as.numeric(median.table[,7])),]
    median.table = data.frame(median.table)
    
    return(median.table)
}
compare_to_null_TajD = compare_to_null(list.of.target.genes, "TajD", null_distribution_TajD)
compare_to_null_nFWH = compare_to_null(ist.of.target.genes, "nFWH", null_distribution_nFWH)
compare_to_null_EW = compare_to_null(list.of.target.genes, "EW", null_distribution_EW)


#4) Complie all the results

ranking_of_median_difference_TajD_simplified = ranking_of_median_difference_TajD[,c(1,10)]
ranking_of_median_difference_nFWH_simplified = ranking_of_median_difference_nFWH[,c(1,10)]
ranking_of_median_difference_EW_simplified = ranking_of_median_difference_EW[,c(1,10)]

compare_to_null_TajD_simplified = compare_to_null_TajD[,c(1,7)]; compare_to_null_TajD_simplified$percentage_ranking = as.numeric(as.character(compare_to_null_TajD_simplified$percentage_ranking))
compare_to_null_nFWH_simplified = compare_to_null_nFWH[,c(1,7)]; compare_to_null_nFWH_simplified$percentage_ranking = as.numeric(as.character(compare_to_null_nFWH_simplified$percentage_ranking))
compare_to_null_EW_simplified = compare_to_null_EW[,c(1,7)]; compare_to_null_EW_simplified$percentage_ranking = as.numeric(as.character(compare_to_null_EW_simplified$percentage_ranking))

for (i in 1:dim(result.target.simplified)[1]){
    gene_id = as.character(result.target.simplified[i,1])
    result.target.simplified[i,21] = as.numeric(ranking_of_median_difference_TajD_simplified[which(ranking_of_median_difference_TajD_simplified$target_id == gene_id),2])
    result.target.simplified[i,22] = as.numeric(compare_to_null_TajD_simplified[which(compare_to_null_TajD_simplified$gene_id == gene_id),2])
    result.target.simplified[i,23] = as.numeric(ranking_of_median_difference_nFWH_simplified[which(ranking_of_median_difference_nFWH_simplified$target_id == gene_id),2])
    result.target.simplified[i,24] = as.numeric(compare_to_null_nFWH_simplified[which(compare_to_null_nFWH_simplified$gene_id == gene_id),2])
    result.target.simplified[i,25] = as.numeric(ranking_of_median_difference_EW_simplified[which(ranking_of_median_difference_EW_simplified$target_id == gene_id),2])
    result.target.simplified[i,26] = as.numeric(compare_to_null_EW_simplified[which(compare_to_null_EW_simplified$gene_id == gene_id),2])
}
colnames(result.target.simplified)[21:26] = c("TajD_rank", "TajD_rank_against_null", "nFWH_rank", "nFWH_rank_against_null", "EW_rank", "EW_rank_against_null")

write.table(result.target.simplified, file = paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/final_data/",indicator, "_DHEW_extracted_compiled_Nov_2017_cleaned_with_function_and_ranking.txt", sep=""), quote=F, row.names = F, col.names = T)



