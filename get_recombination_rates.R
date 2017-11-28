#Getting the recombination rates
#Nov 8th, 2017
#Joo Hyun Im (ji72)


rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/")

coordinates = read.table("Nexus_coordinates_for_all_genes_Nov17.txt", header=T) #695 x 5
coordinates = unique(coordinates) #694 x 5
Dmel_data = read.table("final_data/Dmels_Nov_2017_MK_results_analyzed_final.txt", header=T) #427 genes
#Dmel_data = read.table("DH/Dmel_DHEW_extracted_compiled_Jul_31_2017_cleaned_with_function.txt", header=T, sep="\t") #478 x 12 --old


length(intersect(Dmel_data$gene_id, coordinates$gene.id)) #478
coordinates_new = coordinates[which(coordinates$gene.id %in% Dmel_data$gene_id),]
dim(coordinates_new)

#Print the following way: X:4344444..4345500
command = data.frame(paste(coordinates_new$chr, ":",coordinates_new$coordinate.1,"..",coordinates_new$coordinate.2, sep=""))
write.table(command, file="command_for_getting_recombination_rates_Nov28_2017.txt", quote=F, row.names = F, col.names = F)

#Get the recombination rate data from https://petrov.stanford.edu/cgi-bin/recombination-rates_updateR5.pl?results=0&release=5.36&query_type=batch&Submit=Submit
recom_data = read.table("Nexus_recombination_rates_Nov28_2017.txt", header=T, sep="\t")
coordinates_new_rec = cbind(coordinates_new, recom_data)
write.table(coordinates_new_rec, file="Nexus_recombination_rates_with_gene_ids_Nov28_2017.txt", quote=F, row.names = F, col.names = F)

par(mfrow=c(3,2))
hist(coordinates_new_rec$Comeron.Endpoint.rate, main="Genome-wide", xlab="recombination rate (cM/Mb)"); abline(v=c(2.32), col="red") #genome average
hist(coordinates_new_rec[which(coordinates_new_rec$chr == "2L"),11], main="2L", xlab="recombination rate (cM/Mb)"); abline(v=c(2.39), col="red")
hist(coordinates_new_rec[which(coordinates_new_rec$chr == "2R"),11], main="2R", xlab="recombination rate (cM/Mb)"); abline(v=c(2.66), col="red")
hist(coordinates_new_rec[which(coordinates_new_rec$chr == "3L"),11], main="3L", xlab="recombination rate (cM/Mb)"); abline(v=c(1.79), col="red")
hist(coordinates_new_rec[which(coordinates_new_rec$chr == "3R"),11], main="3R", xlab="recombination rate (cM/Mb)"); abline(v=c(1.96), col="red")
hist(coordinates_new_rec[which(coordinates_new_rec$chr == "X"),11], main="X", xlab="recombination rate (cM/Mb)"); abline(v=c(2.95), col="red")


#Is the focal genes' recombination rate correlative of control genes' recombination rates?
Dmel_data_ordered = Dmel_data[order(Dmel_data$gene_id),]
coordinates_new_rec_ordered = coordinates_new_rec[order(coordinates_new_rec$gene.id), c(1,11)]
Dmel_data_ordered_with_rec = cbind(Dmel_data_ordered, coordinates_new_rec_ordered[,2])
Dmel_data_ordered_with_rec_focal_only = Dmel_data_ordered_with_rec[which(Dmel_data_ordered_with_rec$target_control_status == "target"),]

par=(mfrow=c(1,1))
rec_correlation = data.frame()
for (i in 1:dim(Dmel_data_ordered_with_rec_focal_only)[1]) {
    focal_gene = as.character(Dmel_data_ordered_with_rec_focal_only[i,1])
    focal_gene_rec = as.numeric(Dmel_data_ordered_with_rec_focal_only[i,15])
    control_gene = Dmel_data_ordered_with_rec[grep(focal_gene, Dmel_data_ordered_with_rec$matching_target),]
    control_gene_rec = mean(control_gene[,15])
    rec_correlation[i,1] = focal_gene_rec; rec_correlation[i,2] = control_gene_rec
}
plot(rec_correlation); cor(rec_correlation$V1, rec_correlation$V2) #0.8863
