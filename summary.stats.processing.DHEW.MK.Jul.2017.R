#Process the recent selection summary statistics
#Using DHEW and only considering the results from DH program (TajD, nFWH, EW, HEW, and DHEW)
#July 21, 2017
#Joo Hyun Im (ji72)

########################################################################################################
#A. Get the data
########################################################################################################

#Delete previous data
rm(list=ls(all=TRUE))

#Indicate which species we're working with:
indicator = "Dmel"
#indicator = "Dsim"

#Get the data
get_the_data = function(indicator){
    
    #Read-in the corresponding data
    if (indicator == "Dmel"){
        result.DHEW = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/DH/Dmel_DHEW_extracted_compiled_Jul_31_2017_cleaned.txt", sep="\t", header=T, stringsAsFactors = F) # genes (the genes without enough control genes have not been excluded)

    }
    else{
        result.DHEW = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/DH/Dsim_DHEW_extracted_compiled_Jul_31_2017_cleaned.txt", sep="\t", header=T, stringsAsFactors = F) #477 genes (the genes without enough control genes have not been excluded)
    }
    
    #Add functional categories
    if (indicator == "Dmel"){
        result.label = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/internalization_list_Dmel_Sep_2017.txt", sep="\t", header=T, stringsAsFactors = F) # genes (the genes without enough control genes have not been excluded)
        
    }
    else{
        result.label = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/internalization_list_Dsim_Sep_2017.txt", sep="\t", header=T, stringsAsFactors = F) #477 genes (the genes without enough control genes have not been excluded)
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
write.table(result.total, file="DH/Dmel_DHEW_extracted_compiled_Jul_31_2017_cleaned_with_function.txt", quote=F, row.names = F, col.names = T)

#Check the number of genes in the data
cat("Count the number of genes in the following species: ", indicator, "\n")
result.target = result.total[which(result.total$target.control.status == "target"),]
result.cont = result.total[which(result.total$target.control.status == "control"),]
cat("The number of target genes for ", indicator, " is: ", dim(result.target)[1], "\n")
cat("The number of control genes for ", indicator, " is: ", dim(result.cont)[1], "\n")

#Check if there are control genes for each target gene
result.cont$matching.target = factor(result.cont$matching.target); matching.target.genes = levels(result.cont$matching.target)
dim(result.target)[1]; length(matching.target.genes) #these should be the same


#Split the genes by function

#Autophagy: Induction, Nucleation, Expansion
result.autophagy.target = result.target[which(result.target$type == "Autophagy"),]
result.autophagy.target.induction = result.autophagy.target[which(result.autophagy.target$subtype == "Induction"),]
result.autophagy.target.nucleation = result.autophagy.target[which(result.autophagy.target$subtype == "Nucleation"),]
result.autophagy.target.expansion = result.autophagy.target[which(result.autophagy.target$subtype == "Expansion"),]
#result.autophagy.target.other = result.autophagy.target[which(result.autophagy.target$subtype == "Other"),]
result.autophagy.cont = result.cont[which(result.cont$type == "Autophagy"),]
result.autophagy.cont.induction = result.autophagy.cont[which(result.autophagy.cont$subtype == "Induction"),]
result.autophagy.cont.nucleation = result.autophagy.cont[which(result.autophagy.cont$subtype == "Nucleation"),]
result.autophagy.cont.expansion = result.autophagy.cont[which(result.autophagy.cont$subtype == "Expansion"),]
#result.autophagy.cont.other = result.autophagy.cont[which(result.autophagy.cont$subtype == "Other"),]

#Both
result.both.target = result.target[which(result.target$type == "Both"),]
result.both.cont = result.cont[which(result.cont$type == "Both"),]

#Phagocytosis: Recognition-receptor, Recognition-non-receptor, Internalization, PhagosomeMaturation/FusionWithLysosome, Unknown
result.phago.target = result.target[which(result.target$type == "Phagocytosis"),]
result.phago.target.receptor = result.phago.target[which(result.phago.target$subtype == "Recognition-receptor"),]
result.phago.target.non.receptor = result.phago.target[which(result.phago.target$subtype == "Recognition-non-receptor"),]
result.phago.target.internalization = result.phago.target[which(result.phago.target$subtype == "Internalization"),]
result.phago.target.mat.fus = result.phago.target[which(result.phago.target$subtype == "PhagosomeMaturation/FusionWithLysosome"),]
#result.phago.target.unk = result.phago.target[which(result.phago.target$subtype == "Unknown"),]
result.phago.cont = result.cont[which(result.cont$type == "Phagocytosis"),]
result.phago.cont.receptor = result.phago.cont[which(result.phago.cont$subtype == "Recognition-receptor"),]
result.phago.cont.non.receptor = result.phago.cont[which(result.phago.cont$subtype == "Recognition-non-receptor,"),]
result.phago.cont.internalization = result.phago.cont[which(result.phago.cont$subtype == "Internalization"),]
result.phago.cont.mat.fus = result.phago.cont[which(result.phago.cont$subtype == "PhagosomeMaturation/FusionWithLysosome"),]
#result.phago.cont.unk = result.phago.cont[which(result.phago.cont$subtype == "Unknown"),]

#All internalization:
result.int.target = rbind(result.autophagy.target,result.phago.target,result.both.target) #66 genes
result.int.cont = rbind(result.autophagy.cont,result.phago.cont,result.both.cont) #240 genes

#Recognition:
result.recog.target = result.target[which(result.target$type == "Recognition"),]
result.recog.cont = result.cont[which(result.cont$type == "Recognition"),]

#Signaling:
result.sig.target = result.target[which(result.target$type == "Signaling"),]
result.sig.target.toll = result.sig.target[which(result.sig.target$subtype == "Toll"),]
result.sig.target.imd = result.sig.target[which(result.sig.target$subtype == "IMD"),]
result.sig.cont = result.cont[which(result.cont$type == "Signaling"),]
result.sig.cont.toll = result.sig.cont[which(result.sig.cont$subtype == "Toll"),]
result.sig.cont.imd = result.sig.cont[which(result.sig.cont$subtype == "IMD"),]

#Recognition & Signaling put together:
result.recog.and.sig.target = rbind(result.recog.target, result.sig.target)
result.recog.and.sig.cont = rbind(result.recog.cont, result.sig.cont)


########################################################################################################
#B. Process summary statistics (recent selection)
#Research question: Do autophagy genes or phagocytosis genes as a group show a recent evidence of selection?
#Technical question: Do autophagy genes or phagocytosis genes as a group differ from their respective control genes in terms of summary statistics?
########################################################################################################

#Method 1: Draw distributions of target and control and perform the KS test
#Result shows that the distribution between target and control are not statistically significantly different from each other (9/26/16, 5:10pm)

#With 1:1 controls:
result.cont$matching.target = factor(result.cont$matching.target)
matching.target.genes = levels(result.target$gene_id)
result.cont.one.gene.list = c()
for (i in 1:length(matching.target.genes)){
    gene.name = matching.target.genes[i]
    result.cont.one.gene = result.cont[which(result.cont$matching.target == gene.name),]
    result.cont.one.gene.list = rbind(result.cont.one.gene.list, result.cont.one.gene[1,])
}
result.autophagy.cont = rbind(result.cont.one.gene.list[which(result.cont.one.gene.list$type == "Autophagy"),], result.cont.one.gene.list[which(result.cont.one.gene.list$type == "Both"),]) #added 'Both' category to 'Autophagy (7/23/2017)
result.phago.cont = result.cont.one.gene.list[which(result.cont.one.gene.list$type == "Phagocytosis"),]
result.recog.and.sig.cont = rbind(result.cont.one.gene.list[which(result.cont.one.gene.list$type == "Recognition"),], result.cont.one.gene.list[which(result.cont.one.gene.list$type == "Signaling"),])

#thetaW
median(result.autophagy.target$thetaW); median(result.autophagy.cont$thetaW)
median(result.phago.target$thetaW); median(result.phago.cont$thetaW)
median(result.recog.target$thetaW); median(result.recog.cont$thetaW)
median(result.sig.target$thetaW); median(result.sig.cont$thetaW)
ks.test(result.autophagy.target$thetaW, result.autophagy.cont$thetaW) #Not significant
ks.test(result.phago.target$thetaW,result.phago.cont$thetaW) #Not significant
ks.test(result.recog.target$thetaW, result.recog.cont$thetaW) #Not significant
ks.test(result.sig.target$thetaW, result.sig.cont$thetaW) #Not significant

#TajD
median(result.autophagy.target$TajD); median(result.autophagy.cont$TajD)
median(result.phago.target$TajD); median(result.phago.cont$TajD)
median(result.recog.target$TajD); median(result.recog.cont$TajD)
median(result.sig.target$TajD); median(result.sig.cont$TajD)
ks.test(result.autophagy.target$TajD, result.autophagy.cont$TajD) #Not significant
ks.test(result.phago.target$TajD, result.phago.cont$TajD) #Not significant
ks.test(result.recog.target$TajD, result.recog.cont$TajD) #Not significant
ks.test(result.sig.target$TajD, result.sig.cont$TajD) #Not significant

#FayW
median(result.autophagy.target$nFWH, na.rm = T); median(result.autophagy.cont$nFWH, na.rm = T)
median(result.phago.target$nFWH, na.rm = T); median(result.phago.cont$nFWH, na.rm = T)
median(result.recog.target$nFWH, na.rm = T); median(result.recog.cont$nFWH, na.rm = T)
median(result.sig.target$nFWH, na.rm = T); median(result.sig.cont$nFWH, na.rm = T)
ks.test(result.autophagy.target$nFWH,result.autophagy.cont$nFWH) #Not significant
ks.test(result.phago.target$nFWH,result.phago.cont$nFWH) #Not significant
ks.test(result.recog.target$nFWH,result.recog.cont$nFWH) #Not significant
ks.test(result.sig.target$nFWH,result.sig.cont$nFWH) #Not significant



#Method 2: Create a distribution of [target - median(control)] by assigning each gene in the "cluster" as target. Then overlay where my true target genes are on this distribution. Sometimes the [target - median(control)] value is true (when the true target gene is assigned as target) and other times the value is false (when one of the control genes is assigned as target).

#Restore the control genes
result.autophagy.cont = result.cont[which(result.cont$type == "Autophagy"),]; result.both.cont = result.cont[which(result.cont$type == "Both"),]
result.phago.cont = result.cont[which(result.cont$type == "Phagocytosis"),]; result.recog.and.sig.cont = rbind(result.recog.cont, result.sig.cont)

#Technical question: Are my target genes evenly distributed across the distribution (null) or skewed towards one tail or the other?
distribution_extreme = function(bio_pathway, test_choice, top_or_bottom, percentage){
    cat("We will create a distribution of [target - median(control)] of ", test_choice, " for ", bio_pathway, " genes in ", indicator, " species.", "\n")
    cat("Are my ", bio_pathway, " genes evenly distributed across the distribution (null) or skewed towards ", top_or_bottom, " ", percentage, "%.", "\n")
    
    #Get the right input data
    if (bio_pathway == "Autophagy"){
        #Put gene.id in the matching.target for result.target
        result.target.w.both = rbind(result.autophagy.target,result.both.target) #Dmel
        result.target.w.both$matching.target = result.target.w.both$gene_id
        result.cont.w.both = rbind(result.autophagy.cont,result.both.cont) #Dmel:
        result.for.cluster = rbind(result.target.w.both, result.cont.w.both)
        result.for.cluster$matching.target = factor(result.for.cluster$matching.target)
        list.of.target.genes = levels(result.for.cluster$matching.target)
    } 
    if (bio_pathway == "Phagocytosis"){
        #Put gene.id in the matching.target for result.target
        result.target$matching.target = result.target$gene_id
        result.for.cluster = rbind(result.target, result.phago.cont)
        result.for.cluster$matching.target = factor(result.for.cluster$matching.target)
        list.of.target.genes = levels(result.for.cluster$matching.target)
    }
    
    #For a given cluster, calculate the [target - median(control)] x the number of genes in the cluster and check the distribution of these values
    median.table = data.frame()
    for (n in 1:length(list.of.target.genes)){
        
        if (test_choice == "TajD"){
            test_number = 5 #TajD
        }
        else if (test_choice == "nFWH"){
            test_number = 7 #nFWH
        }
        
        cluster = result.for.cluster[which(result.for.cluster$matching.target == list.of.target.genes[n]),]
        number.of.genes.in.cluster = dim(cluster)[1]
        median.table.cluster = data.frame()
        for (m in 1:number.of.genes.in.cluster){
            target.data = cluster[m,test_number]; median.control.data = median(cluster[-m,test_number])
            diff.median = (target.data - median.control.data) 
            median.table.cluster[m,1] = list.of.target.genes[n]; median.table.cluster[m,2] = cluster[m,1]; median.table.cluster[m,3]= diff.median
            if (m == 1){
                median.table.cluster[m,4] = "Y" #Assigned target status to a true target gene
            }
            else {
                median.table.cluster[m,4] = "N" #Assigned target status to one of the control genes
            }
        }
        median.table = rbind(median.table, median.table.cluster)
    }
    colnames(median.table)[1] = "true.target.gene"; colnames(median.table)[2] = "gene.set.as.target.gene"
    colnames(median.table)[3] = "processed.median.diff"; colnames(median.table)[4] = "target.status.assigned.to.true.target.gene"
    #hist(median.table$processed.median.diff, breaks=10, ylim = c(0,60))
    #shapiro.test(median.table$processed.median.diff) #Normally distributed (pval is significant)
    
    #Rank genes by processed median diff and look into the top and bottom X%
    if (top_or_bottom == "top"){
        median.table.sorted = median.table[order(median.table$processed.median.diff, decreasing=F),] #top is target < median.control (case is more negative)
        top.percentage = median.table.sorted[c(1:as.integer((dim(median.table.sorted)[1]*percentage))),] #Take top X%: target < median.control
    }
    else{
        median.table.sorted = median.table[order(median.table$processed.median.diff, decreasing=T),] #top is target > median.control (case is more positive)
        top.percentage = median.table.sorted[c(1:as.integer((dim(median.table.sorted)[1]*percentage))),] #Take bottom X%: target > median.control
        
    }
    
    #Set expectation and get observation:
    expected.number.of.target.genes = round((dim(median.table.sorted[which(median.table.sorted$target.status.assigned.to.true.target.gene == "Y"),])[1])*percentage, digits=0)
    expected.number.of.cont.genes = round((dim(median.table.sorted[which(median.table.sorted$target.status.assigned.to.true.target.gene == "N"),])[1])*percentage, digits=0)
    observed.number.of.target.genes = dim(top.percentage[which(top.percentage$target.status.assigned.to.true.target.gene == "Y"),])[1]
    observed.number.of.cont.genes = dim(top.percentage[which(top.percentage$target.status.assigned.to.true.target.gene == "N"),])[1]
    
    #Chi-square test
    #sometimes the number of genes in the expected and observed don't add up
    M = matrix(c(expected.number.of.target.genes, expected.number.of.cont.genes, observed.number.of.target.genes, observed.number.of.cont.genes), ncol=2)
    colnames(M) = c("expected", "observed"); rownames(M) = c("target", "control"); 
    
    cat("The fisher exact test results for the following: ", test_choice, " for ", bio_pathway, " genes in ", indicator, " species.", "\n")
    return(fisher.test(M))
}

#Answers:
distribution_extreme("Autophagy", "TajD", "top", 0.1)$p.value #No, autophagy pathway genes are not enriched at the top 10% level
distribution_extreme("Autophagy", "TajD", "top", 0.2)$p.value #No, autophagy pathway genes are not enriched at the top 20% level
distribution_extreme("Autophagy", "TajD", "bottom", 0.1)$p.value #No, autophagy pathway genes are not enriched at the bottom 10% level
distribution_extreme("Autophagy", "TajD", "bottom", 0.2)$p.value #No, autophagy pathway genes are not enriched at the bottom 20% level
distribution_extreme("Phagocytosis", "TajD", "top", 0.1)$p.value #No, phagocytosis pathway genes are not enriched at the top 10% level
distribution_extreme("Phagocytosis", "TajD", "top", 0.2)$p.value #No, phagocytosis pathway genes are not enriched at the top 20% level
distribution_extreme("Phagocytosis", "TajD", "bottom", 0.1)$p.value #No, phagocytosis pathway genes are not enriched at the bottom 10% level
distribution_extreme("Phagocytosis", "TajD", "bottom", 0.2)$p.value #No, phagocytosis pathway genes are not enriched at the bottom 20% level

distribution_extreme("Autophagy", "nFWH", "top", 0.1)$p.value #No, autophagy pathway genes are not enriched at the top 10% level
distribution_extreme("Autophagy", "nFWH", "top", 0.2)$p.value #No, autophagy pathway genes are not enriched at the top 20% level
distribution_extreme("Autophagy", "nFWH", "bottom", 0.1)$p.value #No, autophagy pathway genes are not enriched at the bottom 10% level
distribution_extreme("Autophagy", "nFWH", "bottom", 0.2)$p.value #No, autophagy pathway genes are not enriched at the bottom 20% level
distribution_extreme("Phagocytosis", "nFWH", "top", 0.1)$p.value #No, phagocytosis pathway genes are not enriched at the top 10% level
distribution_extreme("Phagocytosis", "nFWH", "top", 0.2)$p.value #No, phagocytosis pathway genes are not enriched at the top 20% level
distribution_extreme("Phagocytosis", "nFWH", "bottom", 0.1)$p.value #No, phagocytosis pathway genes are not enriched at the bottom 10% level
distribution_extreme("Phagocytosis", "nFWH", "bottom", 0.2)$p.value #No, phagocytosis pathway genes are enriched at the bottom 20% level


#Method 3: For each functional group (autophagy, phagocytosis, both, humoral), calculate the values of target-median(control), compute mean of these values and create a confidence interval. If the CI contains 0, on average  the targets and controls are not statistically significantly different.
library(boot)
create_confidence_interval = function(bio_pathway, test_choice, num_of_rep){
    cat("Create a CI for the ", test_choice, " statistic in ", bio_pathway, " genes.", "\n")
    
    list.of.target.genes = levels(result.total$matching.target) #list of case genes that have control genes (we got the names from matching case column)
    diff.table = data.frame()
    
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
    
    for (n in 1:length(list.of.target.genes)){
        pattern = list.of.target.genes[n]
        data.for.pattern.target = result.total[grep(pattern, result.total[,1]),]
        data.for.pattern.target.test = data.for.pattern.target[,test_number]
        data.for.pattern.target.name = data.for.pattern.target$gene_name
        data.for.pattern.cont = result.total[grep(pattern, result.total[,17]),] #picking the matching target
        data.for.pattern.cont.test.median = median(data.for.pattern.cont[,test_number], na.rm = T)
        type = data.for.pattern.cont$type[1]; subtype = data.for.pattern.cont$subtype[1]
        diff.median = (data.for.pattern.target.test - data.for.pattern.cont.test.median) # difference between target
        diff.table[n,1] = pattern; diff.table[n,2] = data.for.pattern.target.name; diff.table[n,3] = diff.median; diff.table[n,4] = type; diff.table[n,5] = subtype
    } #case minus control
    colnames(diff.table) = c("target.id", "target.name", "diff.median", "type", "subtype") #89 target genes
    diff.table.sorted = diff.table[order(diff.table$diff.median, decreasing=T),]
    diff.table.function = diff.table.sorted[which(diff.table.sorted$type == bio_pathway),]
    cat("Median: ", median(diff.table.function$diff.median, na.rm = T), "\n")
    
    #Create a confidence interval
    #Check the normality assumption
    cat("P-value for checking the normality assumption (Shapiro): ", shapiro.test(c(diff.table.function$diff.median))$p.value, "\n") #If <0.05, NOT NORMALLY DISTRIBUTED
    #Results: Can't use a CI if the data is normally distributed (#http://www.cyclismo.org/tutorial/R/confidence.html)
    
    #Utilize a Bootstrap method in order to calculate a Bootstrap confidence interval
    #Re-sample with replacement B times, for each of these samples calculate the sample mean, and calculate CI.
    #http://stats.stackexchange.com/questions/112829/how-do-i-calculate-confidence-intervals-for-a-non-normal-distribution
    
    #Below is a way to construct a Bootstrap CI, but I have a sample size of <50 (which is not recommended according to Kevin from CSCU)
    #function to obtain the mean
    Bmean = function (data, indices) {
        d = data[indices]; return(mean(d))
    } #goes with library(boot)

    #bootstrapping with X (given) replications
    median_list = na.omit(diff.table.function$diff.median)
    results = boot(data=c(median_list), statistic = Bmean, R=num_of_rep)
    plot(results)
    return(boot.ci(results, type = c("basic"))) #The intervals calculated using the basic bootstrap method.
    #return(boot.ci(results, type = c("norm", "basic", "perc", "bca")))
    #Results: all four CIs for autophagy contain 0, indicating that on average the targets and controls aren't statistically significantly different.    

}

create_confidence_interval("Autophagy", "thetaW", 1000) #pval Sig
create_confidence_interval("Autophagy", "TajD",1000) #pval NS -- normality
create_confidence_interval("Autophagy", "nFWH", 1000) #pval NS -- normality
create_confidence_interval("Autophagy", "EW", 1000)  #pval Sig
create_confidence_interval("Phagocytosis", "thetaW", 1000)
create_confidence_interval("Phagocytosis", "TajD", 1000)
create_confidence_interval("Phagocytosis", "nFWH", 1000)
create_confidence_interval("Phagocytosis", "EW", 1000)

########################################################################################################
#C. Process summary statistics (recent selection)
#Research question: Do sub-categories of autophagy genes or phagocytosis genes show a recent evidence of selection?
########################################################################################################

#Method 1: Simple target-control comparison by t-test (TajD only)
#Recent selection at the functional group level
library(ggplot2); library(grid); library(gridExtra)
library(RColorBrewer)
#display.brewer.all(n=3, type="all", select=NULL, exact.n=TRUE,colorblindFriendly=T) #a list of color combinations that are colorblind-friendly
color_panel = brewer.pal(3, "Dark2") #Call out a panel of 3 colors (minimum 3)

#target + control (median) - autophagy and phagocytosis sub-function
result.int.target = rbind(result.autophagy.target, result.phago.target, result.both.target)
result.int.cont = rbind(result.autophagy.cont, result.phago.cont, result.both.cont)
result.int.target$gene_id = factor(result.int.target$gene_id)
list.of.target.genes = levels(result.int.target$gene_id) #list of int target genes that have control genes (we got the names from matching target column)

#Merge 'Both - Nucleation-PhagosomeMaturation/FusionWithLysosome' into Autophagy 'Other'
# result.int.target$type = levels(revalue(result.int.target$type, c("Both"="Autophagy")))
# result.int.target$subtype = levels(revalue(result.int.target$subtype, c("Nucleation-PhagosomeMaturation/FusionWithLysosome"="Other")))
# result.int.cont$type = levels(revalue(result.int.cont$type, c("Both"="Autophagy")))
# result.int.cont$subtype = levels(revalue(result.int.cont$subtype, c("Nucleation-PhagosomeMaturation/FusionWithLysosome"="Other")))



library(plyr)
median.table = data.frame()
for (n in 1:length(list.of.target.genes)){
    pattern = list.of.target.genes[n]
    data.for.pattern.target = result.total[grep(pattern, result.total[,1]),]
    data.for.pattern.target.TajD = data.for.pattern.target$TajD
    data.for.pattern.cont = result.total[grep(pattern, result.total[,17]),] #matching.target
    #print(dim(data.for.pattern.cont))
    data.for.pattern.cont.TajD.median = median(data.for.pattern.cont$TajD)
    type = data.for.pattern.cont$type[1]; subtype = data.for.pattern.cont$subtype[1]
    median.table[n,1] = pattern; median.table[n,2] = data.for.pattern.target.TajD; median.table[n,3] = data.for.pattern.cont.TajD.median; median.table[n,4] = type; median.table[n,5] = subtype
}
colnames(median.table) = c("target.name", "target.TajD", "cont.TajD.median", "type", "subtype") #Dmel: 72, Dsim: 71 target genes


median.table.rs = reshape(median.table, direction="long", varying = 2:3, idvar = "ids", ids = median.table[,1], timevar="status", v.names="TajD.value", time=c("target","control"))
median.table.rs$subtype <- factor(median.table.rs$subtype, levels = c("Induction", "Nucleation", "Expansion", "Other","Nucleation-PhagosomeMaturation/FusionWithLysosome", "Recognition-receptor", "Recognition-non-receptor", "Internalization", "PhagosomeMaturation/FusionWithLysosome", "Unknown"))
median.table.rs$subtype= revalue(median.table.rs$subtype, c("PhagosomeMaturation/FusionWithLysosome"="Degradation", "Nucleation-PhagosomeMaturation/FusionWithLysosome"="Degradation (both)"))
median.table.rs$status= revalue(median.table.rs$status, c("target"="Focal", "control" = "Control"))
median.table.rs$status <- factor(median.table.rs$status, levels=unique(median.table.rs$status))

#plot with the outline for each point
ggplot(median.table.rs, aes(x=status, y=TajD.value, group=status, color=status))+ scale_color_manual(values = color_panel) + geom_point(aes(color=status), size=12) + geom_point(shape=1, size =12, color = "black") + labs(x = "Functional Category", y="Tajima's D") + ylim(-2.5,1.5) + facet_wrap(~subtype, ncol=5) +theme(axis.text=element_text(size=rel(1.5)), axis.title=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.line=element_line(colour = "black", size=1)) + theme(legend.text=element_text(size=rel(1.5)), legend.title=element_text(size=rel(1.5)), legend.position="top") + theme(strip.text.x = element_text(size = rel(1.5), colour = "black")) + theme_light(base_size = 20)


#t.test: autophagy, both, phagocytosis
name.list = c("Induction", "Nucleation", "Expansion", "Other", "Nucleation-PhagosomeMaturation/FusionWithLysosome","Recognition-receptor", "Recognition-non-receptor","Internalization", "PhagosomeMaturation/FusionWithLysosome", "Unknown") 
t.p.value = data.frame()
for (i in 1:10){
    name = name.list[i]
    median.table.ind = median.table[which(median.table$subtype == name),]; t = t.test(median.table.ind[,2], median.table.ind[,3])
    t.p.value[i,1] = name; t.p.value[i,2] = t$p.value
}
t.p.value
#Answer: Only the Phagosome Maturation/Fusion With Lysosome in phagocytosis genes are different from control genes. After multiple testing correction, 


#Method2: Set up a linear model for all of them (added on 10/14/16) -- comparison between functional groups
#Add the target gene names to matching.target
result.total[1:103,17] = result.total[1:103,1] #for Dmel
#result.total[1:96,17] = result.total[1:96,1] #for Dsim
result.total = result.total[which(result.total$type == "Autophagy" | result.total$type == "Phagocytosis" | result.total$type == "Both"),]
result.total$subtype = droplevels(result.total$subtype)

library(lme4); library(lmerTest); library(lsmeans)
#with interaction:
model = lmer(thetaW ~ target.control.status + subtype + target.control.status:subtype + (1|matching.target), data=result.total); summary(model); anova(model)
model = lmer(TajD ~ target.control.status + subtype + target.control.status:subtype + (1|matching.target), data=result.total); summary(model); anova(model)
model = lmer(nFWH ~ target.control.status + subtype + target.control.status:subtype + (1|matching.target), data=result.total); summary(model); anova(model)
model = lmer(EW ~ target.control.status + subtype + target.control.status:subtype + (1|matching.target), data=result.total); summary(model); anova(model)


#to check for the assumptions of the model:
hist(resid(model)) #expected: centered around 0, normal distribution. looks fine
plot(predict(model),resid(model)) #expected: random scatter around 0, looks okay

lsmeans(model, ~target.control.status) #model means for each subtype
ls1 = lsmeans(model, pairwise~subtype); cld(ls1) #overall comparison of subtypes regardless of target/control status**********
ls2 = lsmeans(model, pairwise~target.control.status|subtype); cld(ls2) #pairwise comparison of target vs control within subtype**********
ls3 = lsmeans(model, pairwise~subtype|target.control.status); cld(ls3) #comparing subtypes within targets (comparing green dots) and within controls (comparing yellow dots)**********

#Re-run in Summer 2017 revealed no statistical significance.
#Previous commmets:Within control genes, there's no sig diff between subtypes. Biologically speaking, this means that these control genes regardless of subtype can be treated the same way. Within target genes, Recognition-receptor genes and PhagosomeMaturation/FusionWithLysosome genes are sig different from each other (but I shouldn't say that Recognition-receptor genes and PhagosomeMaturation/FusionWithLysosome are different from other subtypes). The reason why Recognition-receptor came out as different even if it's ranked in the middle for lsmean may be because it has very low variance.

#taking out the interaction
model2 = lmer(TajD ~ target.control.status + subtype + (1|matching.target), data=result.int.total)
summary(model2)
anova(model2)
#deciding which subtype is different from other subtypes
library(lsmeans); library(multcompView)
lsmeans(model2, ~subtype) #model means for each subtype
ls2 = lsmeans(model2, pairwise~subtype); cld(ls2) #pairwise comparison of subtypes via Tukey test



#Research question: Are there any individual genes that are undergoing selection?
########################################################################################################

#Method 1. Based on TajD or nFWH, rank a target gene compared to the all control genes and see where it falls
#Assumption: The TajD distribution of target and control genes are not statistically significantly different.
#No longer used as of 8/26/2017
rank_target_against_cont_dis = function(result.cont, result.target, statistic){
    if (statistic == "TajD"){
        null_dist = result.cont$TajD
    }
    else if (statistic == "nFWH") {
        null_dist = result.cont$nFWH
    }
    else if (statistic == "EW") {
        null_dist = result.cont$EW
    }
    
    number = as.numeric(length(null_dist))
    sig.target = c(); sig.post.target = c()
    result.target$gene_id = factor(result.target$gene_id)
    list.of.target.genes = levels(result.target$gene_id)
    
    for (i in 1:length(list.of.target.genes)){ #all target genes
        
        if (statistic == "TajD"){
            target_value = result.target[i,5] #TajD
        }
        else if (statistic == "nFWH") {
            target_value = result.target[i,7] #nFWH
        }
        
        all_dist = c(target_value, null_dist); value_label = data.frame(c("target", rep("cont",number)))
        all_dist = cbind(all_dist, value_label); colnames(all_dist) = c("statistic", "status")
        
        ordered_value = all_dist[order(all_dist$statistic),]; rownames(ordered_value) <- seq(length=nrow(ordered_value))
        rank = as.numeric(rownames(ordered_value[which(ordered_value$status == "target"),]))
        
        #top 5% and bottom 5%
        if (rank < ceiling(number*0.05)){
            result = cbind(result.target[i,], rank, as.numeric((rank/number)*100))
            sig.target = rbind(sig.target, result)
        }
        if (rank > (number - ceiling(number*0.05))){ #10% # for positive TajD value ones
            result = cbind(result.target[i,], rank, as.numeric((rank/number)*100))
            sig.post.target = rbind(sig.post.target, result)
        }
    }
    colnames(sig.target)[18] = c(paste("rank.out.of.", number, sep="")); colnames(sig.target)[19] = c("percentage")
    colnames(sig.post.target)[18] = c(paste("rank.out.of.", number, sep="")); colnames(sig.post.target)[19] = c("percentage")
    sig.target.ordered = rbind(sig.target, sig.post.target)
    sig.target.ordered = sig.target.ordered[order(sig.target.ordered$percentage),]
    
    return(sig.target.ordered)

}
ranking_of_target_genes_among_control_genes_TajD = rank_target_against_cont_dis(result.cont, result.target, "TajD")
ranking_of_target_genes_among_control_genes_nFWH = rank_target_against_cont_dis(result.cont, result.target, "nFWH")



#Method 2: Calculate the difference between median of control genes and target gene and rank them.
#list.of.target.genes = levels(result.total$matching.target) #list of case genes that have control genes (we got the names from matching case column)
result.target.internalized = result.target[which(result.target$type %in% c("Autophagy", "Phagocytosis", "Both")),] #only get internalization genes
list.of.target.genes = droplevels(result.target.internalized$gene_id)

rank_the_difference = function(list.of.target.genes, statistic){
    diff.table = data.frame()
    for (n in 1:length(list.of.target.genes)){
        pattern = list.of.target.genes[n]
        data.for.pattern.target = result.total[grep(pattern, result.total[,1]),]
        
        if (statistic == "TajD"){
            data.for.pattern.target.statistic = data.for.pattern.target$TajD
        }
        else if (statistic == "nFWH"){
            data.for.pattern.target.statistic = data.for.pattern.target$nFWH
        }
        else if (statistic == "EW"){
            data.for.pattern.target.statistic = data.for.pattern.target$EW
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
        
    } #target minus control
    
    colnames(diff.table) = c("target_id", "target_name", "target_value", "med_cont_value", "diff_median", "type", "subtype")
    diff.table.sorted = diff.table[order(diff.table$diff_median, decreasing=F),]
    
    return(diff.table.sorted)
}
ranking_of_median_difference_TajD = rank_the_difference(list.of.target.genes, "TajD")
ranking_of_median_difference_nFWH = rank_the_difference(list.of.target.genes, "nFWH")
ranking_of_median_difference_EW = rank_the_difference(list.of.target.genes, "EW")


#Method 3: Based on the p-values (only look at the internalization genes)
result.target.simplified = result.target[which(result.target$type %in% c("Autophagy", "Phagocytosis", "Both")),c(1,2,4,5,6,7,8,9,10,11,12,13,15,16)]
result.target.simplified[,15] = p.adjust(result.target.simplified$TajD_pval, method = "BH")
result.target.simplified[,16] = p.adjust(result.target.simplified$nFWH_pval, method = "BH")
result.target.simplified[,17] = p.adjust(result.target.simplified$EW_pval, method = "BH")
result.target.simplified[,18] = p.adjust(result.target.simplified$DH_pval, method = "BH")
result.target.simplified[,19] = p.adjust(result.target.simplified$HEW_pval, method = "BH")
result.target.simplified[,20] = p.adjust(result.target.simplified$DHEW_pval, method = "BH")
colnames(result.target.simplified)[15:20] = c("TajD_pval_adj", "nFWH_pval_adj", "EW_pval_adj", "DH_pval_adj", "HEW_pval_adj", "DHEW_pval_adj")

result.target.significant.DHEW = result.target.simplified[which(result.target.simplified$DHEW_pval < 0.05),] #22 in Dmel; 9 in Dsim
#Dsim - Worth noting; Vps4, Atg8b, Atg13

result.target.significant.compiled = result.target.simplified[which(result.target.simplified$TajD_pval_adj < 0.05 | result.target.simplified$nFWH_pval_adj < 0.05 | result.target.simplified$EW_pval_adj < 0.05 | result.target.simplified$DH_pval_adj < 0.05 | result.target.simplified$HEW_pval_adj < 0.05 | result.target.simplified$DHEW_pval_adj < 0.05),] #This works for Dmel, but not for Dsim

#For Dsim:
result.target.significant.compiled = result.target.simplified[which(result.target.simplified$TajD_pval < 0.05 | result.target.simplified$nFWH_pval < 0.05 | result.target.simplified$EW_pval < 0.05 | result.target.simplified$DH_pval < 0.05 | result.target.simplified$HEW_pval < 0.05 | result.target.simplified$DHEW_pval < 0.05),] #without multiple testing


#Put together Method 2 and Method 3:
diff_median_table_compiled =c()
for (i in 1:dim(result.target.significant.compiled)[1]){
    gene_id = result.target.significant.compiled[i,1]
    diff_median_TajD = ranking_of_median_difference_TajD[grep(gene_id, ranking_of_median_difference_TajD$target_id),5]
    diff_median_nFWH = ranking_of_median_difference_nFWH[grep(gene_id, ranking_of_median_difference_TajD$target_id),5]
    diff_median_EW = ranking_of_median_difference_EW[grep(gene_id, ranking_of_median_difference_TajD$target_id),5]
    diff_median = cbind(diff_median_TajD, diff_median_nFWH, diff_median_EW)
    diff_median_table_compiled = rbind(diff_median_table_compiled, diff_median)
}
final_sig_cases = cbind(result.target.significant.compiled, diff_median_table_compiled)
#write.table(final_sig_cases, file ="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/individual_genes/Dmel_individual_genes_sig_Aug_26_2017.txt", quote=F, row.names = F, col.names = T)
write.table(final_sig_cases, file ="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/individual_genes/Dsim_individual_genes_sig_but_no_mtc_Aug_26_2017.txt", quote=F, row.names = F, col.names = T)




#Put together the genes that I will evaluate further (8/23/2017)
#Get the following statistics: thetaW, TajD, TajD corrected pval, nFWH, nFWH corrected pval, EW, EW corrected pval, HEW corrected pval
#Dmel
Dmel_evaluate = result.target.simplified[which(result.target.simplified$gene_name %in% c("crq", "YKT6", "CG1599", "Atg8a", "Atg8b")), c(1,2,12,13,3,4,14,6,15,8,16,17,18)]
Dmel_departure_data = cbind(rbind(ranking_of_median_difference_TajD[which(ranking_of_median_difference_TajD$target_name == "crq"),c(3:5)], ranking_of_median_difference_TajD[which(ranking_of_median_difference_TajD$target_name == "CG1599"),c(3:5)], ranking_of_median_difference_TajD[which(ranking_of_median_difference_TajD$target_name == "Atg8b"),c(3:5)], ranking_of_median_difference_TajD[which(ranking_of_median_difference_TajD$target_name == "Atg8a"),c(3:5)], ranking_of_median_difference_TajD[which(ranking_of_median_difference_TajD$target_name == "YKT6"),c(3:5)]), rbind(ranking_of_median_difference_nFWH[which(ranking_of_median_difference_nFWH$target_name == "crq"),c(3:5)], ranking_of_median_difference_nFWH[which(ranking_of_median_difference_nFWH$target_name == "CG1599"),c(3:5)], ranking_of_median_difference_nFWH[which(ranking_of_median_difference_nFWH$target_name == "Atg8b"),c(3:5)], ranking_of_median_difference_nFWH[which(ranking_of_median_difference_nFWH$target_name == "Atg8a"),c(3:5)], ranking_of_median_difference_nFWH[which(ranking_of_median_difference_nFWH$target_name == "YKT6"),c(3:5)]))
colnames(Dmel_departure_data) = c("TajD_target_value", "TajD_med_cont_value", "TajD_diff", "nFWH_target_value", "nFWH_med_cont_value", "nFWH_diff")
Dmel_evaluate_final = cbind(Dmel_evaluate, Dmel_departure_data)
write.table(Dmel_evaluate_final, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/individual_genes/Dmel_individual_genes.txt", quote=F, col.names = T, row.names = F)

#Dsim
Dsim_evaluate = result.target.simplified[which(result.target.simplified$gene_name %in% c("CG1599", "Atg8a", "Atg8b")), c(1,2,12,13,3,4,14,6,15,8,16,17,18)]
Dsim_departure_data = cbind(rbind(ranking_of_median_difference_TajD[which(ranking_of_median_difference_TajD$target_name == "CG1599"),c(3:5)], ranking_of_median_difference_TajD[which(ranking_of_median_difference_TajD$target_name == "Atg8b"),c(3:5)], ranking_of_median_difference_TajD[which(ranking_of_median_difference_TajD$target_name == "Atg8a"),c(3:5)]), rbind(ranking_of_median_difference_nFWH[which(ranking_of_median_difference_nFWH$target_name == "CG1599"),c(3:5)], ranking_of_median_difference_nFWH[which(ranking_of_median_difference_nFWH$target_name == "Atg8b"),c(3:5)],ranking_of_median_difference_nFWH[which(ranking_of_median_difference_nFWH$target_name == "Atg8a"),c(3:5)]))
colnames(Dsim_departure_data) = c("TajD_target_value", "TajD_med_cont_value", "TajD_diff", "nFWH_target_value", "nFWH_med_cont_value", "nFWH_diff")
Dsim_evaluate_final = cbind(Dsim_evaluate, Dsim_departure_data)
write.table(Dsim_evaluate_final, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/individual_genes/Dsim_individual_genes.txt", quote=F, col.names = T, row.names = F)

################################################
#MK TEST
################################################


rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/MK_test_Nov_2016/")

###I. Basic processing
indicator = "Dmel"
#indicator = "Dsim"

basic_processing = function(indicator){
    
    #1. Before reading in the file, I changed the file names listed in the document below to just have the FBgn names.
    #2. Read in the data
    read_data <- function(indicator){
        if (indicator == "Dmel"){
            return(read.table("Dmels_Jul_2017_MK_results_cleaned.txt", header=F, stringsAsFactors = F))
        } else if (indicator == "Dsim"){
            return(read.table("Dsims_Jul_2017_MK_results_cleaned.txt", header=F, stringsAsFactors = F))
        }
    }
    mk_test_results = read_data(indicator)
    colnames(mk_test_results) = c("gene_id","Pn", "Dn", "Ps", "Ds","MKcodons","FETpval")
    number_of_genes = dim(mk_test_results)[1]
    
    
    #3. Add the gene and function names
    read_func_data <- function(indicator){
        if (indicator == "Dmel"){
            return(read.table("internalization_list_Dmel_Sep_2017_MK_only_edited.txt", header=T, sep ='\t', stringsAsFactors = F))
        } else if (indicator == "Dsim"){
            return(read.table("internalization_list_Dsim_Sep_2017_MK_only_edited.txt", header=T, sep='\t', stringsAsFactors = F))
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
    print(ordered_and_selected)
    #results: Rel was the only one with the significant p-value.
    
    #7. Save the data
    write_data <- function(output_data){
        if (indicator == "Dmel"){
            return(write.table(output_data, file="Dmels_Jul_2017_MK_results_analyzed_final.txt", quote=F, col.names = T, row.names = F, sep="\t"))
        } else if (indicator == "Dsim"){
            return(write.table(output_data, file="Dsims_Jul_2017_MK_results_analyzed_final.txt", quote=F, col.names = T, row.names = F, sep="\t"))
        }
    }
    write_data(mk_test_results_final_with_corrected_pval_and_DoS)
    
}
basic_processing(indicator)

###II. Plotting
indicator = "Dmel"
#indicator = "Dsim"
# 
# setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/MK_test_Nov_2016/")
# mk = read.table("internalization_list_Dsim_Sep_2017_MK_only.txt", header=T, sep="\t")
# adj = read.table("../internalization_list_Dsim_Sep_2017.txt", header=T, sep="\t")
# adj_mk =c()
# for (i in 1:dim(mk)[1]){adj_mk = rbind(adj_mk, adj[which(adj$gene_id == as.character(mk[i,1])),])}
# 
# colnames(adj_mk) = c("gene_id", "gene_name", "target_control_status", "type", "subtype", "matching_target")
# write.table(adj_mk, file="internalization_list_Dsim_Sep_2017_MK_only_edited.txt", quote=F, row.names = F, col.names = T, sep = "\t")


#1. Read in data and prepare datasets for plotting:
input_data = read.table(paste(indicator,"s_Jul_2017_MK_results_analyzed_final.txt", sep=""), header=T, sep='\t')
target_only= input_data[which(input_data$target_control_status == "target"),]
target_only = target_only[order(target_only$gene_function),]
target_only_int= target_only[which(target_only$gene_function == "Autophagy" | target_only$gene_function == "Phagocytosis" | target_only$gene_function == "Both"),]
control_only= input_data[which(input_data$target_control_status == "control"),]
control_only = control_only[order(control_only$gene_function),]
control_only_int= control_only[which(control_only$gene_function == "Autophagy" | control_only$gene_function == "Phagocytosis" | control_only$gene_function == "Both"),]

#Only pick the significant cases in the target of internalization genes:
#DoS>0 indicates recurrent positive selection while DoS<0 indicates the presence of purifying selection.
five_percent_sig_pval = 0.05
target_significant_cases_int_only = target_only_int[which(target_only_int$DoS >0 & target_only_int$FETpval < five_percent_sig_pval),]
#Pick the significant cases in the control of all internalization genes:
control_significant_cases_int_only = control_only_int[which(control_only_int$DoS >0 & control_only_int$FETpval < five_percent_sig_pval),]
#Dmel: 0 out of 48 (0%) for target and 11 out of 286 (3.8%) for all controls and 5 out of 206 (2.4%) for int controls
#Dsim: 6 out of 48 (12.5%) for target and 37 out of 286 (12.9%) for all controls and 28 out of 205 (13.7%) for int controls.

#split by function
target_only_int_auto = target_only_int[which(target_only_int$gene_function %in% c('Autophagy', 'Both')),]
target_only_int_phago = target_only_int[which(target_only_int$gene_function %in% c('Phagocytosis')),]
target_significant_cases_int_auto_only = target_only_int_auto[which(target_only_int_auto$DoS >0 & target_only_int_auto$FETpval < five_percent_sig_pval),]
target_only_int_auto[which(target_only_int_auto$DoS <0 & target_only_int_auto$FETpval < five_percent_sig_pval),]
target_significant_cases_int_phago_only = target_only_int_phago[which(target_only_int_phago$DoS >0 & target_only_int_phago$FETpval < five_percent_sig_pval),]

#After FDR correction:
target_significant_cases_adj = target_only_int[which(target_only_int$DoS >0 & target_only_int$corrected_pval < five_percent_sig_pval),]
control_significant_cases_int_only = control_only_int[which(control_only_int$DoS >0 & control_only_int$corrected_pval < five_percent_sig_pval),]


#2. Plot
library(RColorBrewer); library(ggplot2)
list_of_colors_for_function = brewer.pal(5,"Accent")
color_panel = brewer.pal(3, "Dark2") #Call out a panel of 3 colors (minimum 3)
#iii) Plot DoS by boxplot -- target along with its control:
#DoS: In the presence of positive selection, DoS values are positive; values not different from zero are consistent with neutral divergence

#DoS of functions:
input_data$gene_function <- factor(input_data$gene_function, levels =rev(c("Autophagy", "Both", "Phagocytosis", "Recognition", "Signaling")))
input_data$target_control_status <- factor(input_data$target_control_status, levels =c("target","control"))
DoS_plot_all = function(){
    DoS_plot_all = ggplot(aes(x = gene_function, y=DoS), data=input_data) + geom_boxplot(aes(fill= target_control_status), outlier.shape = 1, outlier.color = "black") + coord_flip() + theme(legend.title = element_text(size=15, face="bold"), legend.position="top") + ylab("Direction of Selection (DoS)") + xlab("Function") + geom_hline(yintercept = 0, col="green") + ylim(c(-1,1)) + theme_dark(base_size = 20) + scale_fill_manual(values=color_panel)
    return(DoS_plot_all)
} #The graph is flipped hence the xlab is actually ylab, etc
DoS_plot_all()


#DoS CI
library(boot)
result.total=input_data
create_confidence_interval = function(bio_pathway, test_choice){
    cat("Create a CI for the ", test_choice, " statistic in ", bio_pathway, " genes.", "\n")
    
    list.of.target.genes = levels(result.total$matching_target) #list of case genes that have control genes (we got the names from matching case column)
    diff.table = data.frame()
    
    if(test_choice == "TajD"){
        test_number = 5 #TajD
    }
    else if (test_choice == "nFWH"){
        test_number = 7
    }
    else if (test_choice == "thetaW"){
        test_number = 4 
    }
    else if (test_choice == "DoS"){
        test_number = 14
    }
    
    for (n in 1:length(list.of.target.genes)){
        pattern = list.of.target.genes[n]
        data.for.pattern.target = result.total[grep(pattern, result.total[,1]),]
        data.for.pattern.target.test = data.for.pattern.target[,test_number]
        data.for.pattern.target.name = data.for.pattern.target$gene_name
        data.for.pattern.cont = result.total[grep(pattern, result.total[,12]),] #picking the matching target
        data.for.pattern.cont.test.median = median(data.for.pattern.cont[,test_number], na.rm = T)
        type = data.for.pattern.cont$gene_function[1]; subtype = data.for.pattern.cont$sub_function[1]
        diff.median = (data.for.pattern.target.test - data.for.pattern.cont.test.median) # difference between target
        diff.table[n,1] = pattern; diff.table[n,2] = data.for.pattern.target.name; diff.table[n,3] = diff.median; diff.table[n,4] = type; diff.table[n,5] = subtype
    } #case minus control
    colnames(diff.table) = c("target.id", "target.name", "diff.median", "type", "subtype") #89 target genes
    diff.table.sorted = diff.table[order(diff.table$diff.median, decreasing=T),]
    diff.table.function = diff.table.sorted[which(diff.table.sorted$type == bio_pathway),]
    cat("Median: ", median(diff.table.function$diff.median, na.rm = T), "\n")
    
    #Create a confidence interval
    #Check the normality assumption
    cat("P-value for checking the normality assumption (Shapiro): ", shapiro.test(c(diff.table.function$diff.median))$p.value, "\n") #NS
    #Results: Can't use a CI based on normal distribution, because the data isn't normally distributed (#http://www.cyclismo.org/tutorial/R/confidence.html)
    
    #Utilize a Bootstrap method in order to calculate a Bootstrap confidence interval
    #Re-sample with replacement B times, for each of these samples calculate the sample mean, and calculate CI.
    #http://stats.stackexchange.com/questions/112829/how-do-i-calculate-confidence-intervals-for-a-non-normal-distribution
    
    #Below is a way to construct a Bootstrap CI, but I have a sample size of <50 (which is not recommended according to Kevin from CSCU)
    #function to obtain the mean
    Bmean = function (data, indices) {
        d = data[indices]; return(mean(d))
    } #goes with library(boot)
    
    #bootstrapping with 1000 replications
    median_list = na.omit(diff.table.function$diff.median)
    results = boot(data=c(median_list), statistic = Bmean, R=1000)
    plot(results)
    return(boot.ci(results, type = c("basic"))) #The intervals calculated using the basic bootstrap method.
    #return(boot.ci(results, type = c("norm", "basic", "perc", "bca")))
    #Results: all four CIs for autophagy contain 0, indicating that on average the targets and controls aren't statistically significantly different.    
    
}
create_confidence_interval("Autophagy", "DoS")
create_confidence_interval("Phagocytosis", "DoS")


#DoS REML
library(lme4); library(lmerTest); library(lsmeans)
#Add the target gene names to matching.target
result.total=input_data
result.total = result.total[order(result.total$target_control_status, decreasing = T),]
table(result.total$target_control_status)
#result.total[1:100,12] = result.total[1:100,1] #for Dmel
result.total[1:93,12] = result.total[1:93,1] #for Dsim
result.total = result.total[which(result.total$gene_function == "Autophagy" | result.total$gene_function == "Phagocytosis" | result.total$gene_function == "Both"),]
result.total$sub_function = droplevels(result.total$sub_function)
model = lmer(DoS ~ target_control_status + sub_function + target_control_status:sub_function + (1|matching_target), data=result.total); summary(model); anova(model) #with interaction
ls2 = lsmeans(model, pairwise~target_control_status|sub_function); cld(ls2) #pairwise comparison of target vs control within subtype**********
