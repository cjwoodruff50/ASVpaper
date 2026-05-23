# Rscript to calculate kmer spectra for all ASVs from a dataset, and their 
# associated reads.  Only for execution with many cores (e.g. 56) if there
# are many reads - e.g. 100,000. 
# Also only relevant for mKB and mSerS datasets.
#
# 13 December 2025                                                            [cjw]

args = commandArgs(trailingOnly=TRUE)

#######################################################################################
##########################                               ##############################
##########################         INITIALISATIONS       ##############################
##########################                               ##############################
#######################################################################################
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
packages <- c("stringr","seqinr","tictoc","ShortRead","parallel",
              "stringdist","uwot","rnndescent")    # Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
  #  BiocManager::install(packages[!installed_packages],repos="https://cran.ms.unimelb.edu.au/")
}# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
verbose = FALSE


whichSubMock = 0
if (length(args>0)){ 
  basepath = args[1]
  whichMock = args[2]
  whichSubunit = args[3]
  whichCase = args[4]
  whichSubMock = args[5]
  whichPair = args[6]
  uppErrRate = args[7]
  numcores = args[8]
} else { 
  basepath = "/vast/projects/rrn/ASVtest"
  whichMock = "mockKB"
  whichSubunit = "rrn"
  whichCase = "11"
  whichSubMock = 0
  whichPair = 0
  uppErrRate = "0.01"
  numcores = 2
}


cat("\n\n Code  ASV_reads_kmer_spectra.R starting, with 8 input parameters.\n")
cat(basepath,"\n",whichMock,"\n",whichSubunit,"\n",whichCase,"\n",whichSubMock,"\n",
    whichPair,"\n",uppErrRate,"\n",numcores,"\n\n")

# The error rate is represented by a string with 4 digits to the right of the decimal point. This
# control of representation is necessary for subsequent file naming. 
uppERstr = as.character(uppErrRate)
zz = c2s(sapply(1:6,function(j){out=ifelse(j>nchar(uppERstr),"0",s2c(uppERstr)[j])}))
uppERstr = zz
ERstring = paste("ER","less",uppERstr,sep="_")

if (whichMock=="mockKB"){
  mstr = whichMock 
  dataset = paste("mKB",whichSubunit,paste("C",whichCase,sep=""),sep="_")
  stem1 = paste(mstr,whichSubunit,"Case",whichCase,sep="_")
} else if (whichMock=="mSerF"){
  mstr=whichMock
  dataset = paste(mstr,whichSubunit,paste("C",whichCase,sep=""),sep="_")
  stem1 = paste(mstr,whichSubunit,"Case",whichCase,sep="_")
} else if ("mSerS" == substr(whichMock,1,5)){
  mstr="mSerS"
  dataset = paste(paste(mstr,whichSubMock,sep=""),whichSubunit,paste("C",whichCase,sep=""),sep="_")
  stem1 = paste(whichMock,whichSubunit,"Case",whichCase,sep="_")
} else if (substr(whichMock,1,6) == "mSriSA"){
  mstr="SA"
  dataset = paste(paste(mstr,whichPair,sep=""),whichSubunit,paste("C",whichCase,sep=""),sep="_")
  stem1 = paste(whichMock,whichSubunit,"Case",whichCase,sep="_")
} else if (substr(whichMock,1,6) == "mSriSZ"){
  mstr="SZ"
  dataset = paste(paste(mstr,whichPair,sep=""),whichSubunit,paste("C",whichCase,sep=""),sep="_")
  stem1 = paste(whichMock,whichSubunit,"Case",whichCase,sep="_")
} else {
  mstr ="invalid"
}

cat("stem1:    ",stem1,"    dataset:  ",dataset,"\n")

#######################################################################################
##########################                               ##############################
##########################           FUNCTIONS           ##############################
##########################                               ##############################
#######################################################################################
kmerSpect = function(seq,k){
  # Generate a kmer spectrum for the vector of DNA base characters in seq for 
  # kmers of length k.
  # 26 November 2023                                                [cjw]
  base2b = data.frame(base=c("A","C","G","T"),val=c(0,1,2,3))
  modulo = 4^k
  spectrum = rep(0,4^k)
  mask = 4^k - 1
  j=0
  x=0
  while(j<k){
    j=j+1
    x=4*x;  x=x+base2b[which(base2b$base==seq[j]),"val"]
    #    spectrum[x] = spectrum[x] + 1
  }
  while(j<length(seq)){
    j=j+1
    x=4*x;  x=x+base2b[which(base2b$base==seq[j]),"val"]
    x = x %% modulo
    spectrum[1+x] = spectrum[1+x] + 1
  }
  out = spectrum
}

kmerSpectpar = function(i,readsFastq,kS){
  seq = s2c(as.character(sread(readsFastq[i])))
  out = kmerSpect(seq,kS)
}

#######################################################################################
############################     RUN INITIALISATION      ##############################
#######################################################################################
GRONDpath = file.path(basepath,"GROND")

fastqPath = file.path(basepath,dataset,"fastq")

cat("\n\n    RUNNING   ASV_reads_kmerSpectra.R   \n\n")
cat("rDNA gene: ",whichSubunit,"    Case: ",whichCase, "\n\n")
cat("stem1:    ",stem1,"    dataset:  ",dataset,"\n")

nmock = ifelse(whichMock=="mock10",10,ifelse(whichMock=="mock50",50,ifelse(whichMock=="mockKB",59,7)))

kS = 6      # kmer length used in distance calculations.

#######################################################################################
#######################################################################################
##########################                               ##############################
##########################           MAIN BODY           ##############################
##########################                               ##############################
#######################################################################################
#######################################################################################
cat("\n\nProcessing ",stem1,"\n\n")
# Load   designStrainNames,kbgood,ord.abKB,jord.abund,KB,specR,KBspecR,RA,iSpec,indOp
# inname = paste(dataset,"_characterisation.RData",sep="")
# load(file=file.path(basepath,"RData",inname))


# PART 1: Import text files generated from Julia code call of denoise() followed by 
#         writing of ASV indices to text file.

textpath = file.path(basepath,"text")
RDatapath = file.path(basepath,"RData")
plotpath = file.path(basepath,"plots")
filtfastqpath = file.path(basepath,"fastq")
# Following line only relevant for the Sereika sub-sampling datasets - mstr=="mSerS".
stem2 = "filtered_denoise"
taskSet = c("indices","names","proj","templates")
currentrunID = dataset

# Jindices gives the indices into the filtered fastq file (that is fed to RAD) that
# are associated with each of the ASVs - so Jindices[37] is a string of multiple 
# indices separated by \t associated with ASV 37.  
# Jnames gives the headers of each of the fastq records in the filtered fastq file.
# P is the matrix of 2D UMAP coordinates of the ASVs and reads.
# Templ gives the sequences of each of the ASVs. 
for (whichtask in c(1,2,4)){
  task = taskSet[whichtask]
  if (whichtask==1){
    indname = paste(stem1,ERstring,stem2,task,"out.txt",sep="_")
    Jindices = readLines(con=file.path(textpath,indname))
    cat("Completed reading  Jindices. \n")
  } else if (whichtask==2){
    namesname = paste(stem1,ERstring,stem2,task,"out.txt",sep="_")
    Jnames = readLines(con=file.path(textpath,namesname))
    cat("Completed reading  Jnames.\n")
  } else if (whichtask==3){
    umapname = paste(stem1,ERstring,stem2,task,"out.txt",sep="_")
    P = read.table(file=file.path(textpath,umapname),sep="\t")
    cat("Completed reading  umap projections, P.\n")
  } else {
    templname = paste(stem1,ERstring,stem2,task,"out.txt",sep="_")
    Templ = readLines(con=file.path(textpath,templname))
    cat("Completed reading  templates.\n")
  }
}        #    end    whichtask    loop
nASV = length(Jindices)
nfiltr = length(Jnames)


# Compute kmer spectra for ASVs and all reads and combine in matric kmSpecNormed.
# Note that this should only be done in a slurm run that uses 50+  cores if there
# is a large number of reads.


# Compute kmer spectra for ASVs and all reads and combine in matriX kmSpecNormed.
# Note that this should only be done in a slurm run that uses 50+  cores if there
# is a large number of reads.
# PART 4 subsamples the reads before computing edit distances and computing the 2D UMAP
# calculation.       
Templ.kmer = sapply(1:nASV,function(j){kmerSpect(s2c(Templ[j]),kS)})
# Norm the spectra by dividing by the length of the ASV.  
# We get the length of the ASV by summing the values of its kmer spectrum).
Templ.kmerNormed = sweep(Templ.kmer,MARGIN=2,FUN="/", STATS=colSums(Templ.kmer))   #   Templ.kmer/colSums(Templ.kmer)
# Load the reads fastq files (filtered), extract the sequence data of each read, 
# and compute kmer spectra of these.
if (mstr %in% c("SA","SZ")){
  filtfqName = paste(paste(mstr,whichPair,"_",whichSubunit,"_Case_",whichCase,sep=""),"_",ERstring,"_filtered.fastq",sep="")
} else {
  filtfqName = paste(stem1,ERstring,"filtered.fastq",sep="_")
}
readsFastq = readFastq(file.path(basepath,"fastq",filtfqName))
cat("\n Processing file ",filtfqName, " to generate kmerspectra of",nfiltr, "fastq reads.\n\n")
tic("kmerSpectra computation: ")
ED1 = mclapply(1:nfiltr,kmerSpectpar,readsFastq,kS,mc.cores=numcores)
toc()
rm(readsFastq)
reads.kmer = matrix(unlist(ED1), nrow=4^kS, ncol=nfiltr)
rm(ED1)
hist(colSums(reads.kmer))
reads.kmerNormed = sweep(reads.kmer,MARGIN=2,FUN="/",STATS=colSums(reads.kmer))  
rm(reads.kmer)
nsamp = nfiltr         
kmSpecNormed = cbind(Templ.kmerNormed,reads.kmerNormed)
rm(Templ.kmerNormed,reads.kmerNormed)
outname1 = paste("kmSpecNormed_",stem1,"_",ERstring,".RData",sep="")
# load(file=file.path(basepath,outname1))
save(kmSpecNormed,file=file.path(RDatapath,outname1))


cat(" \n RUN COMPLETED \n")

quit(save = "no", status = 0, runLast = TRUE)

