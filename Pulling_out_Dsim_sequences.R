####################################################
#Pulling out Dsim CDS sequences from WGS
#July 22nd, 2016
#Joo Hyun Im (ji72)
####################################################

#delete any previous input
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/PopGen/Dsimulans_Rogers_2016/")

#1. Convert the FBgn names into CG names
int.genes = read.table("../list.of.case.and.cont.genes.in.Dmel.txt", header=T)
annotation = read.table("../fbgn_annotation_ID_fb_2016_03.tsv", header=F, sep="\t", quote = "", row.names=NULL, stringsAsFactors = F)
colnames(annotation) = c("gene.name", "primary_FBgn", "secondary_FBgn", "primary_CG.name", "secondary_CG.name")
annotation = annotation[,c(1,2,4)]

output.mx <- matrix(NA, ncol=2, nrow=dim(int.genes)[1])
for (i in 1:dim(int.genes)[1]){
    gene.to.compare = as.character(int.genes[i,1])
    checking.against.full.list <- annotation[grep(gene.to.compare, annotation$primary_FBgn),]
    output.mx[i,1] <- as.character(checking.against.full.list[1,1]) #gene.name
    output.mx[i,2] <- as.character(checking.against.full.list[1,3]) #CG.name
}

int.genes.total = cbind(int.genes, output.mx); colnames(int.genes.total) = c("FBgn.number", "gene.name", "CG.name")
write.table(int.genes.total, file ="../list.of.case.and.cont.genes.in.Dmel.FBgn.and.CG.txt", quote=F, row.names = F)
#There were several genes that didn't get the right names (listed as 'NA's). I manually filled them in.


############################################################################################################################################################
#2. Read in GFF file and only get the CDS information
dsim.full.gff = read.table("gff/dsim_update.final_cg_cleaned_up.gff.txt")
dsim.cds.gff <- unique(dsim.full.gff[which(dsim.full.gff$V3== 'CDS'),]) 
int.genes.cg = read.table("information_files/list.of.case.and.cont.genes.FBgn.and.CG.for.Dsim.final.txt", header = T) #I changed 'NA's to right names.
list <- read.table("fly_line_list.txt")

missing.genes = c()
wgs.coordinates = matrix(nrow=dim(int.genes.cg)[1], ncol=7); colnames(wgs.coordinates) = c("FBgn.number", "gene.name","CG.name","chr", "strand","start.coordinate","end.coordinate")
for (i in 1:dim(int.genes.cg)[1]){
    #Look for the CG name of a gene in the GFF file
    gene.of.interest = as.character(paste(int.genes.cg[i,4],"-",sep="")); #8/25/16, I changed [i,3] to [i,4] to reflect the CG names that include the initially missing genes.
    dsim.cds.gff.gene.specific = dsim.cds.gff[grep(gene.of.interest, dsim.cds.gff$V9),]
    
    if (dim(dsim.cds.gff.gene.specific)[1] == 0){ #If the gene is missing in the GFF file,
        missing.genes = rbind(missing.genes, gene.of.interest)
        wgs.coordinates[i,1] = as.character(int.genes.cg[i,1]); wgs.coordinates[i,2] = as.character(int.genes.cg[i,2]); wgs.coordinates[i,3] = as.character(int.genes.cg[i,3]); wgs.coordinates[i,4] = "NA"; wgs.coordinates[i,5] = "NA"; wgs.coordinates[i,6] = "NA"; wgs.coordinates[i,7] = "NA";
        #cat("The coordinates of ", gene.of.interest, " is not available.","\n") #move onto the next gene
    }
    #Out of 510 case and control genes, 103 genes are missing in this GFF file (7/25/2016)
    #The missing ones are all from different chromosomes (some from 3R, others from X, etc.)
    #For 103 missing genes, saved in missing.genes, I used DsimVSDmel.orthos to find the unconventional name for
    #the CG name (ex. 6.g309.t1 is CG4944) and updated the "list.of.case.and.cont.genes.in.Dmel.FBgn.and.CG" to have different names for missing genes.
    #This was finished on 8/14/2016 (check Evernote)
    
    else { #If the gene is in the GFF file,
        
        #Check if there's more than one isoform (added on 8/25/2016)
        dsim.cds.gff.gene.specific$V9 = factor(dsim.cds.gff.gene.specific$V9)
        if (length(levels(dsim.cds.gff.gene.specific$V9))>1){#more than one isoform
            max=c(); chosen.isoform=c()
            for (p in 1:length(levels(dsim.cds.gff.gene.specific$V9))){
                isoform = dsim.cds.gff.gene.specific[which(dsim.cds.gff.gene.specific$V9 == levels(dsim.cds.gff.gene.specific$V9)[p]),]
                max=c(max,isoform[as.numeric(dim(isoform)[1]),5])
            }
            chosen.isoform = levels(dsim.cds.gff.gene.specific$V9)[which.max(max)]
            dsim.cds.gff.gene.specific = dsim.cds.gff.gene.specific[which(dsim.cds.gff.gene.specific$V9 == chosen.isoform),]
        }
        
        if (dsim.cds.gff.gene.specific$V7 == '+'){ #plus strand
            gene.c1 <- as.numeric(dsim.cds.gff.gene.specific[1,4]) -1
            gene.c2 <- as.numeric(dsim.cds.gff.gene.specific[dim(dsim.cds.gff.gene.specific)[1],5]) -1
            
            #Download the sequence from Rogers data using the following coordinates
            #cat("Download the WGS with the following coordinates: ", gene.c1, " and ", gene.c2, "\n")
            wgs.coordinates[i,1] = as.character(int.genes.cg[i,1]); wgs.coordinates[i,2] = as.character(int.genes.cg[i,2]); wgs.coordinates[i,3] = as.character(int.genes.cg[i,3]); wgs.coordinates[i,4] = as.character(dsim.cds.gff.gene.specific[1,1]); wgs.coordinates[i,5] = as.character(dsim.cds.gff.gene.specific[1,7]); wgs.coordinates[i,6] = gene.c1; wgs.coordinates[i,7] = gene.c2
            
            coordinates.mx <- matrix(nrow=dim(dsim.cds.gff.gene.specific)[1], ncol=2)
            for (j in 1:dim(dsim.cds.gff.gene.specific)[1]){
                cds.c1 <- as.numeric(dsim.cds.gff.gene.specific[j,4]); cds.c2 <- as.numeric(dsim.cds.gff.gene.specific[j,5])
                cds.converted.c1 <- cds.c1 - gene.c1 -1; cds.converted.c2 <- cds.c2 - gene.c1 
                if (cds.converted.c1 < 0){ #If there is any negative number for coordinates -> set it to 0.
                    cds.converted.c1 =0
                }
                coordinates.mx[j,1] = cds.converted.c1
                
                #cat("The CDS coordinates are from ", cds.converted.c1, " to ", cds.converted.c2); print("-----")
                if (j == dim(dsim.cds.gff.gene.specific)[1]){
                    coordinates.mx[j,2] <- cds.converted.c2 # in the case of last coordinate
                }
                else{
                    coordinates.mx[j,2] <- cds.converted.c2  
                }
            } #for CDS
            
            colnames(coordinates.mx) <- c("cds_coordinate_1", "cds_coordinate_2")
            coordinates.df <- data.frame(coordinates.mx); coordinates.df.new <- data.frame()
            for (k in 1:dim(coordinates.mx)[1]){
                number.to.compare.1 <- as.numeric(coordinates.mx[k,1])
                number.to.compare.2 <- as.numeric(coordinates.mx[k,2])
                checking <- coordinates.df[which(coordinates.df$cds_coordinate_1 == number.to.compare.1 | coordinates.df$cds_coordinate_2 == number.to.compare.2),]
                #print(checking)
                if (dim(checking)[1]>1){ #When there is more than one isoform
                    min.coordinate.1 <- min(checking[,1])
                    #print(min(checking[,1]))
                    max.coordinate.2 <- max(checking[,2])
                    #print(max(checking[,2]))
                }
                else{
                    min.coordinate.1 <- coordinates.df[k,1]
                    max.coordinate.2 <- coordinates.df[k,2]
                }
                
                coordinates.df.new[k,1] <- min.coordinate.1
                coordinates.df.new[k,2] <- max.coordinate.2
            }
            coordinates.df.mx <- as.matrix(unique(coordinates.df.new))
            
            #print("what is being saved")
            #print(coordinates.df.mx)
            
            ##4.Create a bed file so that I could only get the cds sequences
            #Create a column with the right coordinates
            #If there are multiple CDS regions in the whole gene:
            specific.gene <- matrix(rep(t(coordinates.df.mx),21), ncol=2, byrow=T)
            list.multi <- rep(as.matrix(list), each=dim(coordinates.df.mx)[1]) #or the number of CDS regions
            specific.gene.list <-cbind(list.multi,specific.gene,list.multi)
            write.table(specific.gene.list, file = paste("bed_files/",int.genes.cg[i,1],"_",int.genes.cg[i,2],"_",int.genes.cg[i,3],".bed", sep=""), quote=F, col.names= F, row.names=F, sep = "\t")
        }
        
        else { #minus strand
            gene.c1 <- as.numeric(dsim.cds.gff.gene.specific[1,4]) - 1
            gene.c2 <- as.numeric(dsim.cds.gff.gene.specific[dim(dsim.cds.gff.gene.specific)[1],5]) -1
            #Download the sequence from Nexus suing the following coordinates
            #cat("Download the NEXUS sequence with the following coordinates: ", gene.c1, " and ", gene.c2, "\n")
            wgs.coordinates[i,1] = as.character(int.genes.cg[i,1]); wgs.coordinates[i,2] = as.character(int.genes.cg[i,2]); wgs.coordinates[i,3] = as.character(int.genes.cg[i,3]); wgs.coordinates[i,4] = as.character(dsim.cds.gff.gene.specific[1,1]); wgs.coordinates[i,5] = as.character(dsim.cds.gff.gene.specific[1,7]); wgs.coordinates[i,6] = gene.c1; wgs.coordinates[i,7] = gene.c2
            
            coordinates.mx <- matrix(nrow=dim(dsim.cds.gff.gene.specific)[1], ncol=2)
            for (j in 1:dim(dsim.cds.gff.gene.specific)[1]){
                cds.c1 <- as.numeric(dsim.cds.gff.gene.specific[j,4]); cds.c2 <- as.numeric(dsim.cds.gff.gene.specific[j,5])
                cds.converted.c1 <- cds.c1 - gene.c1 -1; cds.converted.c2 <- cds.c2 - gene.c1 
                if (cds.converted.c1 < 0){ #If there is any negative number for coordinates -> set it to 0.
                    cds.converted.c1 =0
                }
                
                coordinates.mx[j,1] <- cds.converted.c1
                if (j == dim(dsim.cds.gff.gene.specific)[1]){
                    coordinates.mx[j,2] <- cds.converted.c2 # in the case of last coordinate
                }
                else{
                    coordinates.mx[j,2] <- cds.converted.c2 
                }
            } #for CDS
            
            #Only choose the longest isoforms (added and verified on 8/7/2015)
            colnames(coordinates.mx) <- c("cds_coordinate_1", "cds_coordinate_2")
            #print("This is coordinates.mx -1 ")
            #print(coordinates.mx)
            coordinates.df <- data.frame(coordinates.mx)
            coordinates.df.new <- data.frame()
            for (k in 1:dim(coordinates.mx)[1]){
                number.to.compare.1 <- as.numeric(coordinates.mx[k,1])
                number.to.compare.2 <- as.numeric(coordinates.mx[k,2])
                checking <- coordinates.df[which(coordinates.df$cds_coordinate_1 == number.to.compare.1 | coordinates.df$cds_coordinate_2 == number.to.compare.2),]
                #print(checking)
                if (dim(checking)[1]>1){ #When there is more than one isoform
                    min.coordinate.1 <- min(checking[,1])
                    #print(min(checking[,1]))
                    max.coordinate.2 <- max(checking[,2])
                    #print(max(checking[,2]))
                }
                else{
                    min.coordinate.1 <- coordinates.df[k,1]
                    max.coordinate.2 <- coordinates.df[k,2]
                }
                
                coordinates.df.new[k,1] <- min.coordinate.1
                coordinates.df.new[k,2] <- max.coordinate.2
            }
            #print("This is coordinates.df.mx -2 ")
            coordinates.df.mx <- as.matrix(unique(coordinates.df.new))
            #print(coordinates.df.mx)
                        
            #Flip the matrix (only for reverse)
            print("before flipping")
            
            if (dim(coordinates.df.mx)[1] > 1){ #Flip only when there is more than one exon
                rotate <- function(x) t(apply(x, 2, rev))
                coordinates.df.mx <- t(rotate(coordinates.df.mx))
                coordinates.df.mx[dim(coordinates.df.mx)[1],2] <- coordinates.df.mx[dim(coordinates.df.mx)[1],2]
                #print("after flipping and right before creating a bedfile: ")
                #print("This is coordinates.df.mx -3 ")
                #print(coordinates.df.mx)
            }
            else{
                coordinates.df.mx[1,2] <- coordinates.df.mx[1,2]
            }
            
            #print("This is coordinates.df.mx -4 ")
            specific.gene <- matrix(rep(t(coordinates.df.mx),21), ncol=2, byrow=T)
            list.multi <- rep(as.matrix(list), each=dim(coordinates.df.mx)[1]) #or the number of CDS regions
            specific.gene.list <-cbind(list.multi,specific.gene,list.multi)
            #print(specific.gene.list)
            
            write.table(specific.gene.list, file = paste("bed_files/",int.genes.cg[i,1], "_", int.genes.cg[i,2], "_", int.genes.cg[i,3],".bed", sep=""), quote=F, col.names= F, row.names=F, sep = "\t")
        }
    }
    
}
write.table(wgs.coordinates, file="wgs.coordinates.from.Dsim.txt", quote=F, row.names = F)


#3.Create a sh script to pull out Dsim sequences (pseudoCDS) -- added on 8/25/2016
wgs.coordinates = data.frame(wgs.coordinates); wgs.coordinates.na.omit = wgs.coordinates[which(wgs.coordinates$chr != "NA"),] #Get rid of genes whose coordinates cannot be found
sh = c(rep("sh loopGetSeqs.sh",each=as.numeric(dim(wgs.coordinates.na.omit)[1]))); symbol = c(rep(">", each = as.numeric(dim(wgs.coordinates.na.omit)[1])));cg.name.with.fa = paste(wgs.coordinates.na.omit$CG.name,".wgs.fa",sep=""); 
wgs.coordinates.na.omit$chr = as.numeric(wgs.coordinates.na.omit$chr); wgs.coordinates.na.omit$start.coordinate = as.character(wgs.coordinates.na.omit$start.coordinate)
wgs.coordinates.na.omit$end.coordinate = as.character(wgs.coordinates.na.omit$end.coordinate); wgs.coordinates.na.omit$CG.name = as.character(wgs.coordinates.na.omit$CG.name)
script = data.frame(cbind(sh, wgs.coordinates.na.omit$chr, wgs.coordinates.na.omit$start.coordinate, wgs.coordinates.na.omit$end.coordinate, wgs.coordinates.na.omit$CG.name, symbol, cg.name.with.fa))
write.table(script, file="loopGetSeqs_for_all_Dsim_genes.sh", quote=F, row.names = F, col.names = F)



