#Processing the MK test results
#November 9th, 2016 (Updated on 11/14/16 and 11/15/16)
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/MK_test_Nov_2016/")

###I. Basic processing
indicator = "Dmel"
#indicator = "Dsim"
basic_processiong = function(){
    
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
            matched_function = rbind(matched_function, match_function[,c(9:12)])
        }
        return(cbind(mk_test_results$gene_id, matched_gene[,c(2:3)], mk_test_results[,c(2:7)],matched_function))
    }
    mk_test_results_final = add_gene_name_and_function()
    colnames(mk_test_results_final) = c("gene_id","gene_name","gene_name_in_CG","Pn", "Dn", "Ps", "Ds","MKcodons","FETpval", "target_control_status", "gene_function", "sub_function","matching_target")
    
    
    #4.Multiple testing correction
    minus_log_pval = -log(mk_test_results_final$FETpval)
    corrected_pval = p.adjust(mk_test_results_final$FETpval, method = "bonferroni", n = length(mk_test_results_final$FETpval))
    mk_test_results_final_with_corrected_pval = cbind(mk_test_results_final, minus_log_pval, corrected_pval)
    
    
    #5. Calculate DoS (DoS = Dn/(Dn + Ds) − Pn/(Pn + Ps))
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
    
    #7. Save the data
    write_data <- function(output_data){
        if (indicator == "Dmel"){
            return(write.table(output_data, file="Dmels_Nov_2016_MK_results_analyzed_final.txt", quote=F, col.names = T, row.names = F))
        } else if (indicator == "Dsim"){
            return(write.table(output_data, file="Dsims_Nov_2016_MK_results_analyzed_final.txt", quote=F, col.names = T, row.names = F))
        }
    }
    write_data(mk_test_results_final_with_corrected_pval_and_DoS)
    
}
five_percent_sig_pval = 0.05
basic_processiong()

###II. Plotting
indicator = "Dmel"
#indicator = "Dsim"

#1. Read in data and prepare datasets for plotting:
input_data = read.table(paste(indicator,"s_Nov_2016_MK_results_analyzed_final.txt", sep=""), header=T)
target_only= input_data[which(input_data$target_control_status == "target"),]
target_only = target_only[order(target_only$gene_function),]
target_only_int= target_only[which(target_only$gene_function == "Autophagy" | target_only$gene_function == "Phagocytosis" | target_only$gene_function == "Both"),]
control_only= input_data[which(input_data$target_control_status == "control"),]
control_only = control_only[order(control_only$gene_function),]
control_only_int= control_only[which(control_only$gene_function == "Autophagy" | control_only$gene_function == "Phagocytosis" | control_only$gene_function == "Both"),]

#Only pick the significant cases in the target of internalization genes:
target_significant_cases = target_only_int[which(target_only_int$DoS >0 & target_only_int$FETpval < five_percent_sig_pval),]
#Pick the significant cases in the control of all genes (including recognition and signaling): 
control_significant_cases = control_only[which(control_only$DoS >0 & control_only$FETpval < five_percent_sig_pval),]
#Pick the significant cases in the control of all genes (including recognition and signaling): 
control_significant_cases_int_only = control_only_int[which(control_only_int$DoS >0 & control_only_int$FETpval < five_percent_sig_pval),]
#Dmel: 0 out of 48 (0%) for target and 11 out of 286 (3.8%) for all controls and 5 out of 206 (2.4%) for int controls
#Dsim: 6 out of 48 (12.5%) for target and 37 out of 286 (12.9%) for all controls and 28 out of 205 (13.7%) for int controls.



#2. Plot
library(RColorBrewer); library(ggplot2)
list_of_colors_for_function = brewer.pal(5,"Accent")

#i) Plot p-values and draw the significance line:
MK_plot = function(){
    MK_plot_based_on_pval = qplot(x = minus_log_pval, y=gene_name, data=target_only, color=gene_function) + theme(legend.title = element_text(size=15, face="bold"), legend.position="top", legend.text = element_text(size = 15, face = "bold")) + xlab("-log(P value)") + ylab("Gene name") + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) + geom_vline(xintercept = -log(five_percent_sig_pval),col="red") + geom_vline(xintercept = -log(bonnferroni_significant_pval), col="blue")
    MK_plot_based_on_pval + scale_color_manual(values=list_of_colors_for_function, breaks= levels(target_only$gene_function), labels= rev(c("Autophagy", "Both", "Phagocytosis", "Recognition", "Signaling")))
    #Red line denotes the 5% significance level; Blue line denotes Bonferroni correction significance threshold.
    return(MK_plot_based_on_pval)
}
MK_plot()

#ii) Plot DoS by boxplot -- target only:
#DoS: In the presence of positive selection, DoS values are positive; values not different from zero are consistent with neutral divergence

#DoS of functions:
target_only$gene_function <- factor(target_only$gene_function, levels =rev(c("Autophagy", "Both", "Phagocytosis", "Recognition", "Signaling")))
DoS_plot_all = function(){
    DoS_plot_all = qplot(x = gene_function, y=DoS, data=target_only, fill=gene_function) + geom_boxplot() + coord_flip() + theme(legend.title = element_text(size=15, face="bold"), legend.position="right") + ylab("Direction of Selection (DoS)") + xlab("Function") + geom_hline(yintercept = 0, col="dark grey") + ylim(c(-1,1))
    DoS_plot_all + scale_fill_manual(values=list_of_colors_for_function, 
                                  breaks= levels(target_only$gene_function),
                                  labels= rev(c("Autophagy", "Both", "Phagocytosis", "Recognition", "Signaling")))
    return(DoS_plot_all)
} #The graph is flipped hence the xlab is actually ylab, etc
DoS_plot_all()

#DoS of sub-functions:
target_only_int$sub_function <- factor(target_only_int$sub_function, levels = c("Induction","Nucleation","Expansion","Other","Nucleation-PhagosomeMaturation/FusionWithLysosome", "Recognition-receptor", "Recognition-non-receptor", "Internalization", "PhagosomeMaturation/FusionWithLysosome", "Unknown"))
DoS_plot_int = function(){
    DoS_plot_int = qplot(x = sub_function, y=DoS, data=target_only_int, fill=gene_function) + geom_boxplot() + coord_flip() + ylab("Direction of Selection (DoS)") + xlab("Sub-function") + geom_hline(yintercept = 0, col="dark grey") + ylim(c(-1,1))
    DoS_plot_int + scale_fill_manual(values=rev(list_of_colors_for_function[3:5]), 
                                     breaks= levels(droplevels(target_only_int$gene_function)),
                                     labels= c("Autophagy", "Both", "Phagocytosis"))
    return(DoS_plot_int)
} #The graph is flipped hence the xlab is actually ylab, etc
DoS_plot_int()

#iii) Plot DoS by boxplot -- target along with its control:
#DoS: In the presence of positive selection, DoS values are positive; values not different from zero are consistent with neutral divergence

#DoS of functions:
input_data$gene_function <- factor(input_data$gene_function, levels =rev(c("Autophagy", "Both", "Phagocytosis", "Recognition", "Signaling")))
DoS_plot_all = function(){
    DoS_plot_all = ggplot(aes(x = gene_function, y=DoS), data=input_data) + geom_boxplot(aes(color = target_control_status),outlier.shape = NA) + coord_flip() + theme(legend.title = element_text(size=15, face="bold"), legend.position="right") + ylab("Direction of Selection (DoS)") + xlab("Function") + geom_hline(yintercept = 0, col="dark grey") + ylim(c(-1,1))
    return(DoS_plot_all)
} #The graph is flipped hence the xlab is actually ylab, etc
DoS_plot_all()


#DoS of sub-functions:
input_data$sub_function <- factor(input_data$sub_function, levels = c("Induction","Nucleation","Expansion","Other","Nucleation-PhagosomeMaturation/FusionWithLysosome", "Recognition-receptor", "Recognition-non-receptor", "Internalization", "PhagosomeMaturation/FusionWithLysosome", "Unknown"))
DoS_plot_int = function(){
    DoS_plot_int = ggplot(aes(x = sub_function, y=DoS), data=input_data) + geom_boxplot(aes(color = target_control_status),outlier.shape = NA) + coord_flip() + theme(legend.title = element_text(size=15, face="bold"), legend.position="right") + ylab("Direction of Selection (DoS)") + xlab("Function") + geom_hline(yintercept = 0, col="dark grey") + ylim(c(-1,1))
    return(DoS_plot_int)
} #The graph is flipped hence the xlab is actually ylab, etc
DoS_plot_int()


#3. A comparison between target and control in terms of DoS
input_data_matched = input_data[order(input_data$target_control_status, decreasing=T),]
input_data_matched[1:66,13] = input_data_matched[1:66,1] #gets a warning but it works (11/15/2016)

attach(input_data); library(lme4); library(lmerTest); library(lsmeans)
#11/15/16: Results written here are for both Dmel and Dsim

#Model for functional groups:
model = lmer(DoS ~ target_control_status + gene_function + target_control_status:gene_function + (1|matching_target), data=input_data_matched)
summary(model)
anova(model)

#to check for the assumptions of the model:
hist(resid(model)) #expected: centered around 0, normal distribution. looks fine
plot(predict(model),resid(model)) #expected: random scatter around 0, looks okay

lsmeans(model, ~ target_control_status) #model means for each subtype
ls2 = lsmeans(model, pairwise~gene_function); cld(ls2) #overall comparison of subtypes regardless of target/control status
#11/15/16: regardless of target/control status, no gene_function is different from each other.
ls2 = lsmeans(model, pairwise~target_control_status|gene_function); cld(ls2) #pairwise comparison of target vs control within gene function
#11/15/16: No sig diff between tartget and control.
ls3 = lsmeans(model, pairwise~gene_function|target_control_status); cld(ls3) #comparing subtypes within targets and within controls
#11/15/16: Within targets, there's no significant difference. Likewise, there's no significant difference within controls.

#Model for sub-functional groups:
input_data_matched_int_only = input_data_matched[which(input_data_matched$gene_function == "Autophagy" | input_data_matched$gene_function == "Phagocytosis" | input_data_matched$gene_function == "Both"),]

model2 = lmer(DoS ~ target_control_status + sub_function + target_control_status:sub_function + (1|matching_target), data=input_data_matched_int_only)
summary(model2)
anova(model2)

#to check for the assumptions of the model:
hist(resid(model2)) #expected: centered around 0, normal distribution. looks fine
plot(predict(model2),resid(model2)) #expected: random scatter around 0, looks okay

lsmeans(model2, ~ target_control_status) #model means for each subtype
ls2 = lsmeans(model2, pairwise~sub_function); cld(ls2) #overall comparison of subtypes regardless of target/control status
#11/15/16: regardless of target/control status, no gene_function is different from each other.
ls2 = lsmeans(model2, pairwise~target_control_status|sub_function); cld(ls2) #pairwise comparison of target vs control within gene function
#11/15/16: No sig diff between tartget and control.
ls3 = lsmeans(model2, pairwise~sub_function|target_control_status); cld(ls3) #comparing subtypes within targets and within controls
#11/15/16: Within targets, there's no significant difference. Likewise, there's no significant difference within controls.




#
result.int.total = result.total[result.total$type %in% c("Autophagy", "Phagocytosis", "Both"),]
result.int.total$subtype = droplevels(result.int.total$subtype)
