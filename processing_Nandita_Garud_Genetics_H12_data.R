#Processing Dr. Nandita Garud's Genetics 2016 data (on NEXUS-DPGP3 Zambia lines and DGRP lines)
#July 11, 2016, updated on Feb 27, 2018
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE))

setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/comparing_to_existing_data/Nandita_Garud_Genetics_2016_data/")

#Whole chromosome view
par(mfrow = c(4,1))
chr2l.h12 = read.table("Figure2_ZI_Chr2L_above_threshold.txt") #column 1: position, column 2: H12
plot(chr2l.h12[,2]~chr2l.h12[,1], main = c("Zambia 2L"), xlab=c("position"), ylab=c("H12"))
#abline(v=c(449021), col="red") #crq
#abline(v=c(13286940), col="red") #Uvrag
abline(v=c(7987113), col="red") #pes
abline(v=c(15888702), col="red") #TepI
chr2r.h12 = read.table("Chr2R_ZI_RhoAboveThreshold.txt") #column 1: position, column 2: H12
plot(chr2r.h12[,2]~chr2r.h12[,1], main = c("Zambia 2R"), xlab=c("position"), ylab=c("H12"))
#abline(v=c(5730017), col="red") #Vamp7
abline(v=c(20864127), col="red") #emp
chr3l.h12 = read.table("Chr3L_ZI_RhoAboveThreshold.txt") #column 1: position, column 2: H12
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


