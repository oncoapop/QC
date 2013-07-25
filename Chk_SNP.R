##################################################
##  Script to get sequence (which can be SNP masked)    
##  around SNV or indels to check primers already made for 			
##	 Targeted resequencing on the MiSeq 				
##  Aparicio Lab WSOP 2013-001 developed by			
##	 Dr Damian Yap , Research Associate				
##	 dyap@bccrc.ca  Version 2.0 (Jul 2013) 			
##################################################


# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); library('BSgenome.Hsapiens.UCSC.hg19')

library(Biostrings)
library("IRanges")
library("GenomicRanges")
library(Rsamtools)
library('BSgenome.Hsapiens.UCSC.hg19')

# Load the latest available dnSNP for the version of bioconductor if not installed 
len <-  length(available.SNPs())
dbSNP <- available.SNPs()[len]
print("Latest SNP database")
dbSNP

SNP <-   installed.SNPs()
print("Installed SNP database")
dbSNP

# Inject the SNPs into the hg19 ref
SNP_Hsapiens <- injectSNPs(Hsapiens, dbSNP)

#################################################
# Directory structure and file names

basedir="/share/lustre/backup/dyap/Projects"
project="Tumour_Evol"
sample="SA494"

workdir=paste(basedir,project, sep="/")
pydir=paste(workdir,"positions",sep="/")

manfile=paste(sample,"sc.AmpliconManifest",sep="")
pyfile="SA494_pyclone_positions.txt"


######################################################
# input #1 which is the pyclone cluster position file
infile1=paste(pydir, pyfile, sep="/")

# Input #3 which is the AmpliconManfest file
infile2=paste(paste(workdir,sample,sep="/"), manfile, sep="/")

# Output file
outfile=paste(paste(workdir,sample,sep="/"), "SNP_Checked.csv", sep="/")

cluster=as.data.frame(read.table(infile1, header=TRUE, stringsAsFactors = FALSE))
manifest=as.data.frame(read.table(infile2, header=TRUE, skip=5, sep="\t", stringsAsFactors = FALSE))

indf <- data.frame(ID = rep("", nrow(cluster)),
                     Clus = rep(0, nrow(cluster)),
                     stringsAsFactors = FALSE)

for (i in seq(nrow(cluster))) {

                    chr <- strsplit(cluster[i,1], split=":")[[1]][2]
                    pos <-strsplit(cluster[i,1], split=":")[[1]][3]
                    clus <- cluster[i,2]

                    indf$ID[i] <- paste(chr,pos,sep="_")
                    indf$Clus[i] <- clus
                    
                            }

outdf <- data.frame(ID = rep("", nrow(indf)),
                     Clus = rep("", nrow(indf)),
                     AmpSNP = rep("", nrow(indf)),
                     LpriSNP = rep("", nrow(indf)),
                     RpriSNP = rep("", nrow(indf)),
                     Lpriseq = rep("None", nrow(indf)),
                     Rpriseq = rep("None", nrow(indf)),
                     Ampliseq = rep("None", nrow(indf)),
                     SNPLpriseq = rep("NA", nrow(indf)),
                     SNPRpriseq = rep("NA", nrow(indf)),
                     SNPAmpliseq = rep("NA", nrow(indf)),                      
                     stringsAsFactors = FALSE)

for (ri in seq(nrow(indf))) {
  
	        	id <- indf$ID[ri]
  					clus <- indf$Clus[ri]
  					chr <- manifest[manifest[,1] %in% c(id),2]
  					
  					if (length(chr)==0) next
  						else {
						start <-  as.numeric(manifest[manifest[,1] %in% c(id),3])
						end <-  as.numeric(manifest[manifest[,1] %in% c(id),4])
						leftlen <-  as.numeric(manifest[manifest[,1] %in% c(id),5])
						rgtlen <-  as.numeric(manifest[manifest[,1] %in% c(id),6])

# Get the amplicons 1. Reference and 2. SNPmasked
 						ampliseq <- as.character(getSeq(Hsapiens,chr,start,end))
 						snpampliseq <- as.character(getSeq(SNP_Hsapiens,chr,start,end))
 
# Get left and right primers 1. Ref 2. SNPmasked 
						leftend <- start + leftlen
						lpriseq <- as.character(getSeq(Hsapiens,chr,start,leftend))
						snplpriseq <- as.character(getSeq(SNP_Hsapiens,chr,start,leftend))

						rightstart <- end - rgtlen
						rpriseq <- as.character(getSeq(Hsapiens,chr,rightstart,end))
						snprpriseq <- as.character(getSeq(SNP_Hsapiens,chr,rightstart,end))
						
						# Testing to see if the sequence are identical, if they are they do not contain SNPs
						
									if (ampliseq == snpampliseq) ampsnp <- "ok" else { ampsnp <- "SNP" }
									if (lpriseq == snplpriseq) lprisnp <- "ok" else { lprisnp <- "SNP" }
									if (rpriseq == snprpriseq) rprisnp <- "ok" else { rprisnp <- "SNP" }
						
						
# writing the output to the dataframe 
			outdf$ID[ri] <- id
                     	outdf$Clus[ri] <- clus
                     	outdf$AmpSNP[ri] <- ampsnp
                     	outdf$LpriSNP[ri] <- lprisnp
                     	outdf$RpriSNP[ri] <- rprisnp
                     	outdf$Lpriseq[ri] <- lpriseq
                     	outdf$Rpriseq[ri] <- rpriseq
                     	outdf$Ampliseq[ri] <- ampliseq
                     	outdf$SNPLpriseq[ri] <- snplpriseq
                     	outdf$SNPRpriseq[ri] <- snprpriseq
                     	outdf$SNPAmpliseq[ri] <- snpampliseq
							
						}
  }
  
# output file
write.csv(outdf, file = outfile)
