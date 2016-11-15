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
        matched_function = rbind(matched_function, match_function[,c(9:11)])
    }
    return(cbind(mk_test_results$gene_id, matched_gene[,c(2:3)], mk_test_results[,c(2:7)],matched_function))
}
mk_test_results_final = add_gene_name_and_function()
colnames(mk_test_results_final) = c("gene_id","gene_name","gene_name_in_CG","Pn", "Dn", "Ps", "Ds","MKcodons","FETpval", "target_control_status", "function", "sub_function")



#4.Multiple testing correction
minus_log_pval = -log(mk_test_results_final$FETpval)
corrected_pval = p.adjust(mk_test_results_final$FETpval, method = "bonferroni", n = length(mk_test_results_final$FETpval))
mk_test_results_final_with_corrected_pval = cbind(mk_test_results_final, minus_log_pval, corrected_pval)



#5. Calculate DoS (DoS = Dn/(Dn + Ds) − Pn/(Pn + Ps))
five_percent_sig_pval = 0.05
bonnferroni_significant_pval = five_percent_sig_pval/number_of_genes
calculate_DoS = function(data){
    DoS_first = (data$Dn) / (data$Dn + data$Ds); DoS_second = (data$Pn) / (data$Pn + data$Ps)
    DoS = DoS_first - DoS_second
    data_with_DoS = cbind(data, DoS)
    
    return(data_with_DoS[order(data_with_DoS$corrected_pval), ])
}
mk_test_results_final_with_corrected_pval_and_DoS = calculate_DoS(mk_test_results_final_with_corrected_pval)



#6. Select the ones with significant p-values
ordered_and_selected = mk_test_results_final_with_corrected_pval_and_DoS[which(mk_test_results_final_with_corrected_pval_and_DoS$corrected_pval < bonnferroni_significant_pval),]
print(ordered_and_selected)
#results: Rel was the only one with the significant p-value.



#7. Prepare datasets for plotting:
target_only= mk_test_results_final_with_corrected_pval_and_DoS[which(mk_test_results_final_with_corrected_pval_and_DoS$target_control_status == "target"),]
colnames(target_only)[11] = "gene_function"
control_only= mk_test_results_final_with_corrected_pval_and_DoS[which(mk_test_results_final_with_corrected_pval_and_DoS$target_control_status == "control"),]
colnames(control_only)[11] = "gene_function"
target_only = target_only[order(target_only$gene_function),]
target_only_int= rbind(target_only[which(target_only$gene_function == "Autophagy"),], target_only[which(target_only$gene_function == "Phagocytosis"),], target_only[which(target_only$gene_function == "Both"),])
control_only = control_only[order(control_only$gene_function),]
control_only_int= rbind(control_only[which(control_only$gene_function == "Autophagy"),], control_only[which(control_only$gene_function == "Phagocytosis"),], control_only[which(control_only$gene_function == "Both"),])

#Only pick the significant cases:
significant_cases = target_only_int[which(target_only_int$DoS >0 & target_only_int$FETpval < five_percent_sig_pval),]


#8. Plot
library(RColorBrewer); library(ggplot2)
list_of_colors_for_function = brewer.pal(5,"Accent")


#i) Plot p-values and draw the significance line:
target_only$gene_function <- factor(target_only$gene_function, levels =rev(c("Autophagy", "Both", "Phagocytosis", "Recognition", "Signaling")))

MK_plot = qplot(x = minus_log_pval, y=gene_name, data=target_only, color=gene_function) + theme(legend.title = element_text(size=15, face="bold"), legend.position="top", legend.text = element_text(size = 15, face = "bold")) + xlab("-log(P value)") + ylab("Gene name") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + geom_vline(xintercept = -log(five_percent_sig_pval),col="red") + geom_vline(xintercept = -log(bonnferroni_significant_pval), col="blue")
MK_plot + scale_color_manual(values=list_of_colors_for_function, breaks= levels(target_only$gene_function), labels= rev(c("Autophagy", "Both", "Phagocytosis", "Recognition", "Signaling")))
#Red line denotes the 5% significance level; Blue line denotes Bonferroni correction significance threshold.

#ii) Plot DoS by boxplot:

#DoS of functions:
#The graph is flipped hence the xlab is actually ylab, etc


DoS_plot_all = qplot(x = gene_function, y=DoS, data=target_only, fill=gene_function) + geom_boxplot() + coord_flip() + theme(legend.title = element_text(size=15, face="bold"), legend.position="right") + ylab("Direction of Selection (DoS)") + xlab("Function") + geom_hline(yintercept = 0, col="dark grey") + ylim(c(-1,1))
DoS_plot_all + scale_fill_manual(values=list_of_colors_for_function, 
                                  breaks= levels(target_only$gene_function),
                                  labels= rev(c("Autophagy", "Both", "Phagocytosis", "Recognition", "Signaling")))


#DoS of sub-functions:
#The graph is flipped hence the xlab is actually ylab, etc
target_only_int$sub_function <- factor(target_only_int$sub_function, levels = c("Induction","Nucleation","Expansion","Other","Nucleation-PhagosomeMaturation/FusionWithLysosome", "Recognition-receptor", "Recognition-non-receptor", "Internalization", "PhagosomeMaturation/FusionWithLysosome", "Unknown"))

DoS_plot_int = qplot(x = sub_function, y=DoS, data=target_only_int, fill=gene_function) + geom_boxplot() + coord_flip() + ylab("Direction of Selection (DoS)") + xlab("Sub-function") + geom_hline(yintercept = 0, col="dark grey") + ylim(c(-1,1))
DoS_plot_int + scale_fill_manual(values=rev(list_of_colors_for_function[3:5]), 
                                 breaks= levels(droplevels(target_only_int$gene_function)),
                                 labels= c("Autophagy", "Both", "Phagocytosis"))
#In the presence of positive selection, DoS values are positive; values not different from zero are consistent with neutral divergence.





#from: http://stackoverflow.com/questions/21310609/ggplot2-box-whisker-plot-show-95-confidence-intervals-remove-outliers/33639652#33639652
#quantiles_95<- function(x) {
    r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
}
#+ stat_summary(data = target_only, fun.data = quantiles_95, geom="boxplot")



#9. Save the data
write_data <- function(output_data){
    if (indicator == "Dmel"){
        return(write.table(output_data, file="Dmels_Nov_2016_MK_results_analyzed.txt", quote=F, col.names = T, row.names = F))
    } else if (indicator == "Dsim"){
        return(write.table(output_data, file="Dsims_Nov_2016_MK_results_analyzed.txt", quote=F, col.names = T, row.names = F))
    }
}
write_data(mk_test_results_final_with_corrected_pval_and_DoS)

