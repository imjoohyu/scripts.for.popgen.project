#Test if EW test is too liberal
#December 7th, 2017
#Joo Hyun Im (ji72)

#Problem: Lots of focal genes are deviated from control genes (the mean difference between target and median of control is not 0). This is the case in both melanogaster and simulans. "If you do the analysis with a control gene defined as “focal” and contrast it to the other control genes, does this effect go away? Knowing that would help distinguish between a liberal test versus a true effect of these internalization genes." (Brian)

#Delete previous data
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/")

library(plyr)

indicator = "Dmel"
#indicator = "Dsim"

result_total = read.table(paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/final_data/",indicator,"_DHEW_extracted_compiled_Nov_2017_cleaned_with_function.txt", sep=""), header=T, sep='\t') 
colnames(result_total)[14] = c("target_control_status"); colnames(result_total)[17] = c("matching_target")

#Put together autophagy and phagocytosis genes
result_total$type = revalue(result_total$type, c("Autophagy" = "Internalization", "Phagocytosis" = "Internalization", "Both" = "Internalization")) #change [category] into 'Internalization'
result_total$type = factor(result_total$type, levels=unique(result_total$type))

#Scramble targets and controls

one_sample_t_test_scrambled = function(bio_pathway, test_choice){
    #cat("Perform a scrambled 1-sample t-test for the ", test_choice, " statistic in ", bio_pathway, " genes.", "\n")
    
    result_total_target = result_total[which(result_total$target_control_status == "target"),]
    result_total_control = result_total[which(result_total$target_control_status == "control"),]
    result_total_target$matching_target = result_total_target$gene_id
    result_total = rbind(result_total_target, result_total_control)
    list_of_target_genes = c(as.character(result_total_target$gene_id)) #list of case genes that have control genes (we got the names from matching case column)
    
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
    
    for (n in 1:length(list_of_target_genes)){
        cluster = result_total[which(result_total$matching_target == list_of_target_genes[n]),]
        cluster = cluster[order(cluster$target_control_status == "target", decreasing = T),] #Make sure that target is on top
        number.of.genes.in.cluster = dim(cluster)[1]
        
        #Scramble randomly
        random_num = sample(seq(2,number.of.genes.in.cluster),1) #pick something else as a focal gene, excluding a true focal gene
        target_data = cluster[random_num,test_number]; median_control_data = median(cluster[-random_num,test_number])
        diff.median = (target_data - median_control_data) # difference between fake target and control
        type = cluster$type[1]; subtype = cluster$subtype[1]
        diff.table[n,1] = list_of_target_genes[n]; diff.table[n,2] = cluster$gene_name[1]; diff.table[n,3] = diff.median; diff.table[n,4] = type; diff.table[n,5] = subtype
        
    } #case minus control
    
    colnames(diff.table) = c("target.id", "target.name", "diff.median", "type", "subtype")
    diff.table.sorted = diff.table[order(diff.table$diff.median, decreasing=T),]
    diff.table.function = diff.table.sorted[which(diff.table.sorted$type == bio_pathway),]
    #cat("Mean: ", mean(diff.table.function$diff.median, na.rm = T), "\n")
    #cat("Number of genes: ", dim(diff.table.function)[1], "\n")
    
    #Create a confidence interval
    #Check the normality assumption
#     cat("P-value for checking the normality assumption (Shapiro): ", shapiro.test(c(diff.table.function$diff.median))$p.value, "\n")
#     if (shapiro.test(c(diff.table.function$diff.median))$p.value < 0.05){
#         cat("Normality: not normally distributed", "\n") #If <0.05, NOT NORMALLY DISTRIBUTED
#     }
#     else {
#         cat("Normality: normally distributed", "\n")
#     }
    
    t_test_result = t.test(diff.table.function$diff.median, mu=0)
    median_list = na.omit(diff.table.function$diff.median)
    #hist(median_list, main = paste("Hist of ", test_choice, sep = ""))
    
    return(t_test_result[3]) #The intervals calculated using the basic bootstrap method.  #Only give p-value
    
    
}


#Internalization
print("Internalization-EW")
int_EW = c()
for (q in 1:100){
    int_EW = rbind(int_EW, one_sample_t_test_scrambled("Internalization", "EW"))
}

#Signaling
print("Signaling-EW")
sig_EW = c()
for (q in 1:100){
    sig_EW = rbind(sig_EW, one_sample_t_test_scrambled("Signaling", "EW"))
}

#Recognition
print("Recognition-EW")
rec_EW = c()
for (q in 1:100){
    rec_EW = rbind(rec_EW, one_sample_t_test_scrambled("Recognition", "EW"))
}

#Signaling and Recognition together (too few genes otherwise)
result_total$type = revalue(result_total$type, c("Signaling" = "SigRecog", "Recognition" = "SigRecog"))
result_total$type = factor(result_total$type, levels=unique(result_total$type))

print("Signaling/Recognition-EW")
sigrec_EW = c()
for (q in 1:100){
    sigrec_EW = rbind(sigrec_EW, one_sample_t_test_scrambled("SigRecog", "EW"))
}

par(mfrow=c(2,2))
hist(as.numeric(as.character(int_EW)), main = "Internalization", breaks=20)
hist(as.numeric(as.character(sig_EW)), main = "Signaling", breaks=20)
hist(as.numeric(as.character(rec_EW)), main = "Recognition", breaks=20)
hist(as.numeric(as.character(sigrec_EW)), main = "Sig+Rec", breaks=20)


#Added on 2/19/18
#Problem: is the global pattern driven by a few genes that are strongly divergent from the null expectation, or is it a cumulative effect of lots of genes? What happens if you exclude the genes that are individually significant and re-run the analysis? What happens if we remove the top few individual genes with the strongest effect and repeating the analysis? Or are there no genes with individually significant effects?

#Identify genes with the top ~10% of EW pval
result_total_target = result_total[which(result_total$target_control_status == "target"),]
result_total_target_EW_sorted = result_total_target[order(result_total_target$EW_pval),]
result_total_target_EW_sorted_top10per = as.character(result_total_target_EW_sorted[1:10,1])

#Remove them
result_total_top10per_removed = result_total[-which(result_total$gene_id %in% result_total_target_EW_sorted_top10per),]
result_total_top10per_removed = result_total_top10per_removed[-which(result_total_top10per_removed$matching_target %in% result_total_target_EW_sorted_top10per),]
result_total_top10per_removed = droplevels(result_total_top10per_removed)

#Re-run the analysis
one_sample_t_test = function(bio_pathway, test_choice){
    cat("Perform a 1-sample t-test for the ", test_choice, " statistic in ", bio_pathway, " genes.", "\n")
    
    list.of.target.genes = levels(result_total_top10per_removed$matching_target) #list of case genes that have control genes (we got the names from matching case column)
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
        list.of.target.genes = levels(result_total_top10per_removed$matching_target)
        matching_target_number =12
    }
    
    for (n in 1:length(list.of.target.genes)){
        pattern = list.of.target.genes[n]
        data.for.pattern.target = result_total_top10per_removed[grep(pattern, result_total_top10per_removed[,1]),]
        data.for.pattern.target.test = data.for.pattern.target[,test_number]
        data.for.pattern.target.name = data.for.pattern.target$gene_name
        data.for.pattern.cont = result_total_top10per_removed[grep(pattern, result_total_top10per_removed[,matching_target_number]),] #picking the matching target
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
    #hist(median_list, main = paste("Hist of ", test_choice, sep = ""))
    
    return(t_test_result) #The intervals calculated using the basic bootstrap method.
    
    
}
one_sample_t_test("Internalization", "EW")


