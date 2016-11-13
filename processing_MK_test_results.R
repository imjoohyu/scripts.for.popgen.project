#Processing the MK test results
#November 9th, 2016
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/MK_test_Nov_2016/")

indicator = "Dmel"
#indicator = "Dsim"

#1. Before reading in the file, I changed the file names listed in the document below to just have the FBgn names.
#2. Read in the data
read_data <- function(indicator){
    if (indicator == "Dmel"){
        return(read.table("Dmels_Nov_2016_MK_results.txt", header=F))
    } else if (indicator == "Dsim"){
        return(read.table("Dsims_Nov_2016_MK_results.txt", header=F))
    }
}
mk_test_results = read_data(indicator)
colnames(mk_test_results) = c("gene_id","Pn", "Dn", "Ps", "Ds","MKcodons","FETpval")
number_of_genes = dim(mk_test_results)[1]

#3. Add the gene and function names
gene_list = read.table("../Dsimulans_Rogers_2016/comparison.between.Dmel.and.Dsim/list.of.case.and.cont.genes.for.Dsim.txt", head =T)
function_list = read.table("../Dsimulans_Rogers_2016/Dsim.summary.statistics.110416.w.selected.control.genes.txt", head =T)
add_gene_name_and_function = function(){
    matched_gene =c(); matched_function =c()
    for (i in 1:dim(mk_test_results)[1]){
        match_gene = gene_list[which(gene_list$FBgn.number == as.character(mk_test_results[i,1])),]
        match_function = function_list[which(function_list[,1]== as.character(mk_test_results[i,1])),]
        matched_gene = rbind(matched_gene, match_gene[,c(1:3)])
        matched_function = rbind(matched_function, match_function[,c(10:11)])
    }
    return(cbind(mk_test_results$gene_id, matched_gene[,c(2:3)], mk_test_results[,c(2:7)],matched_function))
}
mk_test_results_final = add_gene_name_and_function()
colnames(mk_test_results_final) = c("gene_id","gene_name","gene_name_in_CG","Pn", "Dn", "Ps", "Ds","MKcodons","FETpval", "function", "sub_function")

#4.Multiple testing correction
corrected_pval = p.adjust(mk_test_results_final$FETpval, method = "bonferroni", n = length(mk_test_results_final$FETpval))
mk_test_results_final_with_corrected_pval = cbind(mk_test_results_final, corrected_pval)

#5. Calculate DoS (DoS = Dn/(Dn + Ds) − Pn/(Pn + Ps))
significant_pval = 0.05/number_of_genes
calculate_DoS = function(data){
    DoS_first = (data$Dn) / (data$Dn + data$Ds); DoS_second = (data$Pn) / (data$Pn + data$Ps)
    DoS = DoS_first - DoS_second
    data_with_DoS = cbind(data, DoS)
    
    return(data_with_DoS[order(data_with_DoS$corrected_pval), ])
}
mk_test_results_final_with_corrected_pval_and_DoS = calculate_DoS(mk_test_results_final_with_corrected_pval)

#6. Select the ones with significant p-values
ordered_and_selected = mk_test_results_final_with_corrected_pval_and_DoS[which(mk_test_results_final_with_corrected_pval_and_DoS$corrected_pval < significant_pval),]
print(ordered_and_selected)
#results: Rel was the only one with the significant p-value.

#7. Save the data
write_data <- function(output_data){
    if (indicator == "Dmel"){
        return(write.table(output_data, file="Dmels_Nov_2016_MK_results_analyzed.txt", quote=F, col.names = T, row.names = F))
    } else if (indicator == "Dsim"){
        return(write.table(output_data, file="Dsims_Nov_2016_MK_results_analyzed.txt", quote=F, col.names = T, row.names = F))
    }
}
write_data(ordered)

