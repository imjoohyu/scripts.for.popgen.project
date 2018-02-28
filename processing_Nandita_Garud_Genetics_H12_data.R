#Processing Dr. Nandita Garud's Genetics 2016 data (on NEXUS-DPGP3 Zambia lines and DGRP lines)
#July 11, 2016, updated on Feb 27, 2018
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE))

#No data on ChrX
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/comparing_to_existing_data/Nandita_Garud_Genetics_2016_data/")


#Method 1) Check if the coordinates overlap with peaks (Feb 27, 2018)

#A. Pull coordinates of focal genes (start or mid)
target_genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Manuscript/PopGen/figures_tables/Feb2018/Dmel_target_genes_final_Feb2018.txt", head=T)
coordinates = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/Dmelanogaster_Nexus/Dmel_Nexus_rerun_after_fixing_isoform_issue_Sept_23_2016/Nexus_coordinates_for_all_internalization_genes_for_Dmel_Sept_23_2016.txt", head=T)

corresponding_coordinates = c()
corresponding_chr = c()
for (i in 1:dim(target_genes)[1]){
    gene_id = as.character(target_genes[i,1])
    #gene_coordinate = as.numeric(coordinates[which(coordinates$gene.id == gene_id),2]) #start coordinate
    gene_coordinate = 0.5*(as.numeric(coordinates[which(coordinates$gene.id == gene_id),3]) - as.numeric(coordinates[which(coordinates$gene.id == gene_id),2])) + as.numeric(coordinates[which(coordinates$gene.id == gene_id),2]) #coordinate of the middle of the gene 
    
    chr = as.character(coordinates[which(coordinates$gene.id == gene_id),4])
    corresponding_coordinates = rbind(corresponding_coordinates, gene_coordinate)
    corresponding_chr = rbind(corresponding_chr, chr)
}

target_genes_w_coordinates = cbind(target_genes, corresponding_coordinates, corresponding_chr)


#B. Compare the coordinates of list of defined peaks
Zambia_peaks = read.table("TableS4_Top_25_peaks_Zambia.txt", header = T)

overlap_with_defined_peaks = function(chr){
    Zambia_peaks_chr = Zambia_peaks[which(Zambia_peaks$Chromosome == paste("Chr",chr,sep="")),]
    target_genes_chr = target_genes_w_coordinates[which(target_genes_w_coordinates$corresponding_chr == chr),]
    
    save_overlap = c()
    for (m in 1:dim(target_genes_chr)[1]){
        target_gene_coor = as.numeric(target_genes_chr[m,4])
        for (r in 1:dim(Zambia_peaks_chr)[1]){
            if (Zambia_peaks_chr[r,2] <target_gene_coor & Zambia_peaks_chr[r,3] > target_gene_coor){
                print("Overlapped gene: "); print(target_genes_chr[m,])
                cat("Overlapped peak:"); print(Zambia_peaks_chr[r,])
                save_overlap = rbind(save_overlap, target_genes_chr[m,])
            }
        }
    }

    chr_H12 = read.table(paste("Figure2_ZI_Chr", chr, "_above_threshold.txt", sep=""))
    plot(chr_H12[,2]~chr_H12[,1], main = c(paste("Zambia_", chr, sep="")), xlab=c("position"), ylab=c("H12"))
    abline(v=c(save_overlap$corresponding_coordinates), col="red")
    }

par(mfrow = c(4,1))
overlap_with_defined_peaks("2L")
overlap_with_defined_peaks("2R")
overlap_with_defined_peaks("3L")
overlap_with_defined_peaks("3R")



#
#C. A more automatic way to do this: Check coordinates against H12 values -- FAIL

get_coordinates = function(chr){
    chr_H12 = read.table(paste("Figure2_ZI_Chr", chr, "_above_threshold.txt", sep="")) #first column: location, second column: H12
    target_genes_chr = target_genes_w_coordinates[which(target_genes_w_coordinates$corresponding_chr == chr),]
    plot(chr_H12[,2]~chr_H12[,1], main = c(paste("Zambia_", chr, sep="")), xlab=c("position"), ylab=c("H12"))
    
    highest_H12 = c()
    for (j in 1:dim(target_genes_chr)[1]){
        rough_coor_range_1 = round(as.numeric(target_genes_chr[j,4]), -3)
        
        if (rough_coor_range_1 < max(chr_H12$V1) && rough_coor_range_1 > min(chr_H12$V1)){
            rough_coor_range_2 = rough_coor_range_1 + 300000
            chr_H12_selected = chr_H12[which(chr_H12$V1 > rough_coor_range_1 & chr_H12$V1 < rough_coor_range_2),]
            highest_H12_selected = max(chr_H12_selected$V2)
        }
        else {
            highest_H12_selected = NA #coordinate is outside the range of data
        }
        highest_H12 = rbind(highest_H12, highest_H12_selected)
    }
    
    target_genes_chr = cbind(target_genes_chr, highest_H12)
    
    return(target_genes_chr)
    
}
check_peaks =function(target_genes_chr, rough_peak_size){
    if (dim(target_genes_chr[which(target_genes_chr$highest > rough_peak_size),])[1] > 0){
        print("One of the genes is on/close to the probable peak")
    }
    else{
        print("No gene is on/close to the probable peak")
    }
    return(target_genes_chr[which(target_genes_chr$highest > rough_peak_size),])
}
par(mfrow = c(4,1))

coordinates_2L = get_coordinates("2L")
coordinates_2L_peaked = check_peaks(coordinates_2L, 0.035) #0.04+ (manually entered based on the plot)
abline(v=c(coordinates_2L_peaked$corresponding_coordinates), col="red")

coordinates_2R = get_coordinates("2R")
coordinates_2R_peaked =check_peaks(coordinates_2R, 0.045) #0.045+ (based on the plot)
abline(v=c(coordinates_2R_peaked$corresponding_coordinates), col="red")

coordinates_3L = get_coordinates("3L")
coordinates_3L_peaked =check_peaks(coordinates_3L, 0.035) #0.035+ (based on the plot)
abline(v=c(coordinates_3L_peaked$corresponding_coordinates), col="red")
#gets -Inf for park but it's so far out on the side that it doesn't matter

coordinates_3R = get_coordinates("3R")
coordinates_3R_peaked =check_peaks(coordinates_3R, 0.04) #0.035+ (based on the plot)
abline(v=c(coordinates_3R_peaked$corresponding_coordinates), col="red")


#Method 2) Compare the coordinate visually (July 11, 2016)

#Whole chromosome view
par(mfrow = c(4,1))
chr2l.h12 = read.table("Figure2_ZI_Chr2L_above_threshold.txt") #column 1: position, column 2: H12
plot(chr2l.h12[,2]~chr2l.h12[,1], main = c("Zambia 2L"), xlab=c("position"), ylab=c("H12"))
#abline(v=c(449021), col="red") #crq
#abline(v=c(13286940), col="red") #Uvrag
abline(v=c(7987113), col="red") #pes
abline(v=c(15888702), col="red") #TepI
chr2r.h12 = read.table("Figure2_ZI_Chr2R_above_threshold.txt") #column 1: position, column 2: H12
plot(chr2r.h12[,2]~chr2r.h12[,1], main = c("Zambia 2R"), xlab=c("position"), ylab=c("H12"))
#abline(v=c(5730017), col="red") #Vamp7
abline(v=c(20864127), col="red") #emp
chr3l.h12 = read.table("Figure2_ZI_Chr3L_above_threshold.txt") #column 1: position, column 2: H12
plot(chr3r.h12[,2]~chr3r.h12[,1], main = c("Zambia 3R"), xlab=c("position"), ylab=c("H12"))
abline(v=c(7426992), col="red") #Rac2
#abline(v=c(136670), col="red") #Vps24
#abline(v=c(4818982), col="red") #CG11975
chr3r.h12 = read.table("Figure2_ZI_Chr3R_above_threshold.txt") #column 1: position, column 2: H12
plot(chr3r.h12[,2]~chr3r.h12[,1], main = c("Zambia 3R"), xlab=c("position"), ylab=c("H12"))
abline(v=c(4818982), col="red") #CG11975
abline(v=c(24942113), col="red") #Atg14
abline(v=c(25660826), col="red") #Atg16


#Zoomed-in chromosome view
par(mfrow = c(1,2))
chr2l.h12.zoom = chr2l.h12[which(chr2l.h12[,1] > 10000000 & chr2l.h12[,1] < 15000000), ]
plot(chr2l.h12.zoom[,2]~chr2l.h12.zoom[,1], main = c("Zambia 2L"), xlab=c("position"), ylab=c("H12"))
abline(v=c(449021), col="red") #crq
abline(v=c(13286940), col="red") #Uvrag
chr3r.h12.zoom = chr3r.h12[which(chr3r.h12[,1] < 7000000), ]
plot(chr3r.h12.zoom[,2]~chr3r.h12.zoom[,1], main = c("Zambia 3R"), xlab=c("position"), ylab=c("H12"))
abline(v=c(4175856), col="red") #Atg13
abline(v=c(136670), col="red") #Vps24


