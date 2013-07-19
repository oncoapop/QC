# This script QCs the no of reads that can be assigned from the given samplesheet in a run

#############################
RunID="130603_M00897_0033_000000000-A49AR"
Alignment="Alignment"
##############################
Aligndir=paste(paste("/Volumes/Monco/MiSeq Run Files/MiSeq Analysis Files",RunID,sep="/"), "Data/Intensities/BaseCalls", sep="/")
dir=paste(Aligndir,Alignment,sep="/")
setwd(dir)

##############################
outdir="/Users/dyap/Documents/USB_STORAGE/QC"

#Change the skip to make sure the rows are read in correctly (skipping previous section)
counts <- read.table(file="DemultiplexSummaryF1L1.txt", stringsAsFactors=FALSE, sep="\t", skip=34)
samples <- read.csv(file="SampleSheetUsed.csv", stringsAsFactors = FALSE, skip=32, header= TRUE)

#########################################################################
##     DO NOT CHANGE ANYTHING BELOW THIS LINE - WSOP QC PIPELINE      ##

# This reads the run name from the samplesheet and removes space from the name
samplesheet <- read.csv(file="SampleSheetUsed.csv", stringsAsFactors = FALSE)
expt <- samplesheet[4,2]
fname <- paste("Total_Reads_QC", gsub("\\s","", expt), sep="_")

# The names of the outputfiles are according to the Expt name in the sample sheet
pdffile=paste(paste(outdir,fname,sep="/"), "pdf",sep=".")
csvfile=paste(paste(outdir,fname,sep="/"), "csv",sep=".")

# Processing files for comparison
names(counts)[1]<-"Index"
names(counts)[2]<-"IndexRC"
names(counts)[3]<-"Count"


outdf <- data.frame(Name = rep("Unidentified", nrow(counts)),
                     Barcode = rep("", nrow(counts)),
   		   Reads = rep(0, nrow(counts)),
                     stringsAsFactors = FALSE)
                     
# Matches the samples vs reads by barcode 
for (ri in seq(nrow(counts))) {
  tmp <- samples[grep(counts$Index[ri], samples$Index),]

  name <- tmp$Sample_Name
  bc <- counts$Index[ri]
  gp <- tmp$Group
  des <- tmp$Description
  cnt <- as.numeric(counts$Count[ri])
  
if (length(name) >1 ) {outdf$Name[ri] <- "Filtered"} else if (length(name) !=0 ) {outdf$Name[ri] <- name}
if (length(bc) !=0 ) {  outdf$Barcode[ri] <- bc}
  outdf$Reads[ri] <- cnt
 
  }

pdf(pdffile, width=6, height=6)

barplot(outdf$Reads, xlab="Sample", ylab="No of Reads", main="Total Reads per Sample", names.arg=outdf$Name, las=3, cex.names=0.5)

dev.off()

		  									
write.csv(outdf, csvfile)

