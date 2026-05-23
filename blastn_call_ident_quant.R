# Rscript to construct a blastn call given a multi-sequence fasta file of queries and
# a blastn reference database, and then to process blastn text output to extract
# identification and quantification of microbiome strains.
# This is derived from blastn_call_analysis_v03.R with internal date of 23 October 2025.
# This code is written to use the simple directory structure of /vast/projects/rrn/ASVcode,
# and also to be as minimalist as is consistent with the output requirements.
#
# It is written to only process the Srinivas et al. datasets, the Sereika et al. extracted
# datasets, and datasets for the mKB mock microbiome based on King et al.'s Gut Feeling
# Knowledge Base and simulated nanopore reads from Wick's badread code.
#
# source("/vast/projects/rrn/RscriptsArchive/current/blastn_call_ident_quant.R")
#
#
# 23 April 2026                                                   [cjw]

args = commandArgs(trailingOnly=TRUE)
#######################################################################################
##########################                               ##############################
##########################         INITIALISATIONS       ##############################
##########################                               ##############################
#######################################################################################
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
packages <- c("stringr","seqinr","tictoc","readxl","ShortRead","stringdist",
              "uwot","rnndescent","parallel")
# c("stringr","seqinr","tictoc","readxl","ShortRead", "DECIPHER","bfsl","dplyr",
#              "nnet", "MASS","stringdist","compositions","scales","uwot","rnndescent")    # Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Now set values of parameters specifying the microbiome being analysed and of the analysis.
# whichMock can take the following forms:- "mockKB", "mSriSA1","mSriSZ3","mSerF","mSerS2", 
# where the digit, if it occurs, can be 1,2 3 or 4.
if (length(args>0)){
  basepath = args[1]
  whichMock = args[2]
  whichSubunit = args[3]
  whichCase = args[4]
  whichSubMock = args[5]
  whichPair = args[6]
  uppErrRate = args[7]
  maxopnum = args[8]
  nmock = args[9]
  numcores = args[10]
} else {
  basepath = "/vast/projects/rrn/ASVtest"
  whichMock = "mockKB"
  whichSubunit = "rrn"
  whichCase = "11"
  whichSubMock = 0
  whichPair = 0
  uppErrRate = 0.01
  maxopnum = 291  # 261 for mock50, 48 for mock10   
  nmock = 59     #  50 for mock50, 10 for mock10
  numcores = 4
}
cat("\n\n Code blastn_call_ident_quant.R starting, with 10 input parameters.\n")
cat(basepath,"\n",whichMock,"\n",whichSubunit,"\n",whichCase,"\n",whichSubMock,"\n",
             whichPair,"\n",uppErrRate,"\n",maxopnum,"\n",nmock,"\n",numcores,"\n\n")

# The error rate is represented by a string with 4 digits to the right of the decimal point. This
# control of representation is necessary for subsequent file naming. 
uppERstr = as.character(uppErrRate)
zz = c2s(sapply(1:6,function(j){out=ifelse(j>nchar(uppERstr),"0",s2c(uppERstr)[j])}))
uppERstr = zz
ERstring = paste("ER","less",uppERstr,sep="_")
verbose = FALSE

#######################################################################################
##########################                               ##############################
##########################           FUNCTIONS           ##############################
##########################                               ##############################
#######################################################################################

splitNamesSxy = function(id){ # id = S23[[k]]$ID[whichRank]
  t1 = unlist(strsplit(id,split="_"))
  if (length(t1)==1){ # Need to split on space
    t1 = unlist(strsplit(id,split=" "))
  }
  # Remove [ ] or ' 'from those genus names that have them.
  tc = s2c(t1[1])
  if (tc[1] %in% c("[","'")){
    tc = tc[2:(length(tc)-1)]
    t1[1] = c2s(tc)
  }
  genus = t1[1]
  species = paste(t1[1:2],sep=" ", collapse=" ")
  strain = paste(t1[3:(length(t1)-1)],sep=" ", collapse=" ")
  strainfull = paste(t1[1:(length(t1)-1)],sep=" ", collapse=" ")
  operon = t1[length(t1)]
  operonfull = paste(t1,sep=" ", collapse=" ")
  out = c(genus,species,strainfull,operonfull)
}


generate_RA = function(nmock,M){
  # Generate relative abundances forming a geometric series of nmock elements ranging from
  # 1 to M.
  # 29 July 2025                                                                  [cjw]
  r = exp(log(M)/(nmock-1)) - 1
  RA = sapply(1:nmock,function(j){out=M/((1+r)^(j-1))})
  print(RA)
  RA = RA/sum(RA)
  out = RA
}

which_specR = function(specName){
  out=which(sapply(1:length(specR),function(j){specR[[j]]$species.ordered==specName}))
}

which_specRfullC = function(specName){
  out=which(sapply(1:length(specRfullC),function(j){specRfullC[[j]]$species.ordered==specName}))
}

OptoStrain = function(opName){
  t1 = unlist(strsplit(opName,split=" "))
  out = paste(t1[1:(length(t1)-1)],sep=" ",collapse=" ")
}

OptoSpecies = function(opName){
  t1 = unlist(strsplit(opName,split=" "))
  out = paste(t1[1:2],sep=" ",collapse=" ")
}

OptoGenus = function(opName){
  t1 = unlist(strsplit(opName,split=" "))
  out = t1[1]
}

create_specRfullC = function(fastaName){
  # Construct key object, specRfullC.  
  # Prior quantities calculation:-
  #      igs = which(sapply(1:length(opDB.strains),function(j){nchar(opDB.strains[j])>1}))
  #      nigs = length(igs)
  #      specswithstrain = opDB.species[igs]
  #      uniqSpecwithStrain = unique(specswithstrain)
  #      iord1 = order(uniqSpecwithStrain)
  F1 = read.fasta(file=file.path("/vast/projects/rrn/GROND/output/fasta",fastaName),seqtype="DNA",
                  as.string=TRUE,forceDNAtolower = FALSE,whole.header=TRUE)
  opDB.species = sapply(1:length(F1),function(j){attributes(F1[j])$name})
  # For species names need to trim off the index.
  speciesNames = sapply(1:length(opDB.species),function(j){
                            t1 = unlist(strsplit(opDB.species[j],split="_"))
                            out = paste(t1[1:2],sep="_",collapse="_")
  })
  nSpec = length(opDB.species)
  uniqSpeciesNames = unique(speciesNames);    nuniqSpec = length(uniqSpeciesNames)
  specRfullC = vector("list",nuniqSpec)
  countmiss=0
  index.missed = rep(0,100)
  for (j in 1:nuniqSpec){  
    indDB = which(speciesNames==uniqSpeciesNames[j])
    if (length(indDB)>0){
      istr = which(sapply(1:length(indDB),function(j){nchar(speciesNames[indDB[j]]) >1}))
      strainNames = opDB.species[indDB[istr]]
      specRfullC[[j]] = list(idb=indDB[istr], species.ordered=uniqSpeciesNames[j], specOps=strainNames)
      #  print(specR[[j]])
    } else {
      countmiss = countmiss + 1
      index.missed[countmiss] = j
      cat("No member of opDB.species is identical to uniqSpecwithStrain ",j," which is ",uniqSpeciesNames[j],"\n")
    }
  }
  out = specRfullC
}


linKBtospecR = function(grondPath){
  # Only relevant for mockKB that has 59 strains.
  # Creates key objects and mappings for 
  #    1. Identifying which KB strains have entries in specR;
  #    2. Ordering KB strains by relative abundance;
  #    3. Mapping abundance-ordered KB strains to their corresponding specR entries.
  #    4. Determining strain operons for any KB strain in specR.
  # 11 August 2025 
  # Load crapOpNames, opDB.species, opDB.strains, OpNames, ops_per_strain_each_spec, specR
  inname = "specR_etal_grondDB.RData"
  load(file.path(grondPath,inname))
  innameK = "KBspecR_matrix.RData"
  load(file.path(grondPath,innameK))
  kbgood = which(KBspecR[,"opHi"]>0);  nkbgd = length(kbgood)
  # for (j in 1:nkbgd){print(specR[[KBspecR[kbgood[j],1]]]$strainOps[KBspecR[kbgood[j],2]:KBspecR[kbgood[j],3]])}
  
  # Get KB and KB abundance data to guide rank-ordering on KB abundance data of the strains being used 
  # in mock microbiomes.
  gfkb.name = "KGFDB.xlsx"  #  "King CH etal SuppMat S4 Table GutFeeling.KB PLoSOne.0206484.s008.xlsx"
  KB = read_excel(file.path(grondPath,gfkb.name))
  abundKB.name = "KGFDB_Abundance_Tables.xlsx"
  AbundKB = read_excel(file.path(grondPath,abundKB.name))
  ord.abKB = order(AbundKB$Average,decreasing=TRUE)
  # Relation between KB and AbundKB is derived with the following:-
  # Define jord.abund as an index on kbgood such that KB[kbgood[jord.abund]],] gives
  # the KB rows that have strains ordered in descending relative abundances - that is,
  # KB[kbgood[jord.abund[j,1]],1:3] corresponds to AbundKB[ord.abKB[jord.abund[j,2]],1c(1,3,4)]
  jord.abund = matrix(0,nrow=59,ncol=2)
  jord.abund[1,] = c(33,2);    jord.abund[2,] = c(14,8);    jord.abund[3,] = c(34,13)
  jord.abund[4,] = c(22,14);   jord.abund[5,] = c(52,20);   jord.abund[6,] = c(58,21)
  jord.abund[7,] = c(65,24);   jord.abund[8,] = c(24,26);   jord.abund[9,] = c(54,33)
  jord.abund[10,] = c(71,42);  jord.abund[11,] = c(9,44);  jord.abund[12,] = c(28,48)
  jord.abund[13,] = c(85,49);  jord.abund[14,] = c(67,59);  jord.abund[15,] = c(35,60)
  jord.abund[16,] = c(66,61);  jord.abund[17,] = c(25,62);  jord.abund[18,] = c(57,63)
  jord.abund[19,] = c(11,64);  jord.abund[20,] = c(47,65);  jord.abund[21,] = c(73,67)
  jord.abund[22,] = c(60,68);  jord.abund[23,] = c(55,70);  jord.abund[24,] = c(88,72)
  jord.abund[25,] = c(30,73);  jord.abund[26,] = c(59,74);  jord.abund[27,] = c(36,75)
  jord.abund[28,] = c(72,78);  jord.abund[29,] = c(26,79);  jord.abund[30,] = c(10,80)
  jord.abund[31,] = c(18,81);  jord.abund[32,] = c(39,84);  jord.abund[33,] = c(37,85)   
  jord.abund[34,] = c(50,95);  jord.abund[35,] = c(62,98);  jord.abund[36,] = c(79,99)
  jord.abund[37,] = c(46,100); jord.abund[38,] = c(68,106); jord.abund[39,] = c(86,108)
  jord.abund[40,] = c(45,110); jord.abund[41,] = c(61,111); jord.abund[42,] = c(16,113)
  jord.abund[43,] = c(75,115); jord.abund[44,] = c(27,116); jord.abund[45,] = c(31,119)
  jord.abund[46,] = c(64,120); jord.abund[47,] = c(32,121); jord.abund[48,] = c(80,125)
  jord.abund[49,] = c(17,128); jord.abund[50,] = c(51,129); jord.abund[51,] = c(29,130)    
  jord.abund[52,] = c(56,135); jord.abund[53,] = c(84,137); jord.abund[54,] = c(70,139)
  jord.abund[55,] = c(77,144); jord.abund[56,] = c(83,148); jord.abund[57,] = c(89,149)
  jord.abund[58,] = c(63,150); jord.abund[59,] = c(23,155)
  
  ra1 = AbundKB[ord.abKB[jord.abund[,2]],16]
  RAkbg = ra1/sum(ra1)
  out = list(kbgood,ord.abKB,jord.abund,KB,specR,KBspecR,RAkbg)
}

ASVOpCount = function(opInd){
  # Given a vector of operon indices,  returns a table giving the count of
  # each unique index.
  # 24 June2025                                                          [cjw]
  uniq = unique(opInd);   nuniq = length(uniq)
  opCounts = sapply(1:nuniq,function(j){length(which(opInd==uniq[j]))})
  T = data.frame(OpIndex=uniq,opCount=opCounts)
  out = T
}

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

editDistAlgn = function(seq1, seq2){
  # Computes the Levenshtein distance between base sequences seq1, seq2
  s1 = seq1;   L1 = nchar(s1)
  s2 = seq2;   L2 = nchar(s2)
  if ( (nchar(s1)>1100) && (nchar(s2)>1100)){
    algnf = pairwiseAlignment(s1,s2,type="overlap",substitutionMatrix=mat,
                              gapOpening=gapOpen,gapExtension=gapExtend)
    algnr = pairwiseAlignment(s1,reverseComplement(DNAString(s2)),type="overlap",substitutionMatrix=mat,
                              gapOpening=gapOpen,gapExtension=gapExtend)
    if (score(algnf)>score(algnr)){algn=algnf} else {algn=algnr}
    alscore = score(algn)
    PID = pid(algn)
    Lev = stringDist(c(as.character(pattern(algn))[1],as.character(subject(algn))))
  } else {
    alscore = -999;  Lev = -1;  PID = -999
  }
  out = c(alscore, PID, Lev)
}  

plot_quality_histo = function(inputName, plotName,ThreshLow,ThreshHi){
  # Take the pre-filtered or filtered fastq data from julia code for RAD denoising and generate 
  # distributional data on the reads. This is relevant to both real and simulated reads. 
  # First, extract header data from the fastq file fed to RAD and save as a headers_....txt file.
  filtered = length(grep("filtered",inputName))==1
  headersName = paste(paste("headers",stem1,ERstring,ifelse(filtered,"filtered","pre-filtered"),sep="_"),".txt",sep="")
  if (filtered){
    argstrh1 = paste(" '^@seq' ",file.path(basepath,"fastq",inputName),
                     " > ",file.path(basepath,"text",headersName),sep="")
    system2("grep", args=argstrh1) 
    
    Hf = read.delim(file.path(basepath,"text",headersName),sep="|",header=FALSE,col.names=c("seqID","ErrRate"))
    nfiltr = length(Hf$ErrRate)
    errorRate = as.numeric(sapply(1:nfiltr, function(j){unlist(strsplit(Hf$ErrRate[j],split="="))[2]}))
    iret = sapply(1:nfiltr,function(j){as.numeric(substr(Hf$seqID[j],start=5,stop=nchar(Hf$seqID[j])))})
    nr = length(Hf)
    errorRate = as.numeric(sapply(1:nfiltr, function(j){unlist(strsplit(Hf$ErrRate[j],split="="))[2]}))
    iret = sapply(1:nfiltr,function(j){as.numeric(substr(Hf$seqID[j],start=5,stop=nchar(Hf$seqID[j])))})
    
    Hf.df = data.frame(Index=Hf$seqID, ErrorRate=errorRate)
    pdf(file=file.path(basepath,"plots",plotName), paper="a4")
    hist(errorRate,xlab="Estimated Error Rate", ylab="Frequency",
         main=paste("Estimated ErrorRate of Reads of filtered ", dataset,sep=" "))
    dev.off()
    
  } else{
    argstr2 = paste(" '^@' ",file.path(basepath,"fastq",inputName)," > ",file.path(basepath,"text","temp_headers_prefilt.txt"),sep="")
    system2("grep", args = argstr2)
    argstr3 = paste(" 'strand' ",file.path(basepath,"text","temp_headers_prefilt.txt")," > ",file.path(basepath,"text",headersName),sep="")
    system2("grep", args = argstr3)
    Hf = (read.delim(file=file.path(basepath,"text",headersName),sep=" ", header=FALSE))
    nr = length(Hf$V1)
    # Now parse each header element to extract the ONTID OperonName Strand Length and PID.
    
    ONTID = rep("",nr);   OpName = rep("",nr);   Strand = rep("",nr);
    Length = rep(0.0,nr);   PID = rep(0.0,nr)
    for (jr in 1:nr){
      ONTID[jr] = Hf[jr,1]    
      t2 = unlist(strsplit(Hf[jr,2],split=","));    OpName[jr] = t2[1];  
      st = substr(t2[2],start=1,stop=1);   Strand[jr] = ifelse(st=="-",-1,1) 
      Length[jr] = as.numeric(substr(Hf[jr,3], start=8, stop=nchar(Hf[jr,3])));      
      pd = unlist(strsplit(Hf[jr,5],split="[=]"))[2];   PID[jr] = as.numeric(substr(pd,start=1, stop=nchar(pd)-1))
    }
    T.df = data.frame(ID=ONTID, operon=OpName, strand=Strand, readLength=Length, fragLength=Length,PID=PID)
    ifilt1 = which(sapply(1:nr,function(j){out=(T.df$PID[j]>99.00) && (T.df$readLength[j]>ThreshLow && T.df$readLength[j]<ThreshHi)})) 
    jlengthFilt = which(sapply(1:nr,function(j){T.df$readLength[j]>ThreshLow && T.df$readLength[j]<ThreshHi}))
    pdf(file=file.path(basepath,"plots",plotName), paper="a4")
    hist(100-T.df$PID[jlengthFilt],xlab="100-PID", ylab="Frequency", breaks = seq(from=0.0,to=1.1*max(100-T.df$PID[jlengthFilt]), by = 0.05),
         main=paste("Percent Bases In Alignment Not Identical; length, not quality, filtered " , dataset,sep=" "))
    dev.off()
    
  } 
  if(filtered){
    out = iret
  } else {
    out = list(ifilt1 = ifilt1, Tfdf = T.df)
  }
}

editDistabs = function(sp1, sp2, k){
  # Computes the approximate Levenshtein distance between sequences
  # having length-nmormed kmer spectra sp1, sp2.
  #    sp1, sp2     the kmer spectrum for kmers of length  k  of sequences 1 and 2
  #    Lseq1, Lseq2 the length in bases of sequences 1 and 2.
  #  29 June 2025                                              [cjw]
  tot = sum(sapply(1:length(sp1),function(j){out=abs(sp1[j] - sp2[j])}))
  out = tot/(2*k)
}

EDkmerPar = function(i,allSpec,kS){
  # allSpec is a matrix of 4^kS rows, each column being a length-normed kmer spectrum
  # of a read for k-mer length kS.
  # Returns a vector of ncol(allSpecs) edit distances.
  # 29 June 2025                                               [cjw]
  sp1 = allSpec[,i]
  out = sapply(1:ncol(allSpec),function(j){editDistabs(sp1,allSpec[,j],kS)})
}
#######################################################################################
############################     RUN INITIALISATION      ##############################
#######################################################################################
# IMPORTANT: Must ensure the system PATH parameter includes path to the executable of blastn
#     export PATH=$PATH:/vast/projects/rrn/GROND/ncbi-blast-2.16.0+/bin
#  and that BLASTDB is also defined    
#     export BLASTDB=/vast/projects/rrn/GROND/blastdb1
system2("export",args="PATH=$PATH:/vast/projects/rrn/GROND/ncbi-blast-2.16.0+/bin")
system2("export",args="BLASTDB=/vast/projects/rrn/GROND/blastdb1")
# Alignment quality thresholds (simple)
mmthresh16 = 10
gapthresh16 = 5
mmthreshITS = 7
gapthreshITS = 3
mmthresh23 = 17
gapthresh23 = 8
mmthreshrrn = 30
gapthreshrrn = 12


if (whichMock=="mockKB"){
  mstr = whichMock 
  dataset = paste("mKB",whichSubunit,paste("C",whichCase,sep=""),sep="_")
  stem1 = paste(whichMock,whichSubunit,"Case",whichCase,sep="_")
} else if (whichMock=="mSerF"){
  mstr="mSerF"
  dataset = paste(mstr,whichSubunit,paste("C",whichCase,sep=""),sep="_")
  stem1 = paste(whichMock,whichSubunit,"Case",whichCase,sep="_")
} else if ("mSerS" == substr(whichMock,1,5)){
  mstr="mSerS"
  dataset = paste(paste(mstr,whichSubMock,sep=""),whichSubunit,paste("C",whichCase,sep=""),sep="_")
  stem1 = paste(whichMock,whichSubunit,"Case",whichCase,sep="_")
} else if (substr(whichMock,1,6) == "mSriSA"){
  mstr="SA"
  dataset = paste(paste(mstr,whichPair,sep=""),whichSubunit,paste("C",whichCase,sep=""),sep="_")
  stem1 = paste(paste(mstr,whichPair,sep=""),whichSubunit,"Case",whichCase,sep="_")
} else if (substr(whichMock,1,6) == "mSriSZ"){
  mstr="SZ"
  dataset = paste(paste(mstr,whichPair,sep=""),whichSubunit,paste("C",whichCase,sep=""),sep="_")
  stem1 = paste(paste(mstr,whichPair,sep=""),whichSubunit,"Case",whichCase,sep="_")
} else {
  mstr ="invalid"
}

cat("stem1:    ",stem1,"    dataset:  ",dataset,"\n")

grondPath = file.path(basepath,"GROND")
blastDBpath = file.path(grondPath,"blastdb2")
if (whichMock %in% c("mock10","mock50","mockKB","mSerF","mSerS")){
  blastDBname = "rrnstrainDB_grond_refseq207.fasta"
  blastDBname.short = "rrnstrainGrondRefSeq"
  blastDBpath = file.path(grondPath,"blastdb1")
} else { # Covering for the Srinivas datasets.
  blastDBname = "rrnSpeciesDB_grond_refseq207full_C.fasta"
  blastDBname.short = "rrnSpeciesGrondRefSeq"
  blastDBpath = file.path(grondPath,"blastdb2")
}
RADFastaPath = file.path(basepath,"fasta")
outRDatapath = file.path(basepath,"RData")
plotpath = file.path(basepath,"plots")
textpath = file.path(basepath,"text")

ASVfastaName = paste(stem1,"_",ERstring,"_filtered.fasta",sep="")

# Load crapOpNames, opDB.species, opDB.strains, OpNames, ops_per_strain_each_spec, specR
inname = "specR_etal_grondDB.RData"
load(file.path(grondPath,inname))

# Set various booleans that determine what parts of the code below are executed.
doAlign = FALSE    # Set to FALSE if alignment has been done previously as the output from 
# alignment will have been generated.  Also the initial processing that parses this
# output does not need to be repeated.

computeKmerSpectra = FALSE   # only need to be TRUE for first run with this dataset.

fullReadset = FALSE    # If TRUE the full set  of ASVs + reads is used to compute the 2D 
#                       UMAP projection using uwot - see lines ~1441, ~1630.
# The setting of the Boolean, projJulia, controls which method is used.  Setting
# projJulia == TRUE requires that the Julia code providing the 2D UMAP coordinates has been run prior to the 
# plotting code below.
projJulia = ifelse(whichMock == "mockKB",FALSE,TRUE)
# If  prefiltplot  is TRUE plots are generated of read quality data for reads in the primary mock microbiome libary 
# whose alignments were of acceptable length.  This can take some time for large libraries, and is generally not used.
prefiltplot = FALSE 

#######################################################################################
#######################################################################################
##########################                               ##############################
##########################           MAIN BODY           ##############################
##########################                               ##############################
#######################################################################################
#######################################################################################
# PART 0: Sets up the designed mock microbiome variables for the chosen tasks.  This
#         includes the designed relative abundances, the number of operons per strain, 
#         and indOp, the matrix indexing the first and last operon for each strain.
#         NOTE: The number of operons per strain has been named num_ops_per_strains.
if (whichMock %in% c("mock50","mockKB")){
  nmock = ifelse(whichMock=="mock50",50,59)
  # DKB below gives  list(kbgood,ord.abKB,jord.abund,KB,specR,RAkbg)
  DKB = linKBtospecR(grondPath)
  kbgood = DKB[[1]];              nkbgd = length(kbgood)
  ord.abKB = DKB[[2]] 
  jord.abund = DKB[[3]];          njab = nrow(jord.abund)
  KB = DKB[[4]]
  specR = DKB[[5]]
  KBspecR = DKB[[6]]
  RAkbg = DKB[[7]]
  
  designStrainNames = sapply(1:njab,function(j){
    specRinds3 = KBspecR[kbgood[jord.abund[j,1]],]
    t0 = specR[[specRinds3[1]]]$strainOps[specRinds3[2]];  
    t1 = unlist(strsplit(t0,split=" "))
    out = paste(t1[1:(length(t1)-1)],sep=" ",collapse=" ")})
  ind.strainOps = KBspecR[kbgood[jord.abund[,1]],"opLow"]
  indupp.strainOps = KBspecR[kbgood[jord.abund[,1]],"opHi"]
  num_ops_per_strains = sapply(1:length(ind.strainOps),function(j){indupp.strainOps[j] - ind.strainOps[j]+1})
  maxopnum = sum(num_ops_per_strains)
  designStrainNames = designStrainNames[1:nmock]
  
  designSpeciesNames = sapply(1:nmock,function(j){
    specRinds3 = KBspecR[kbgood[jord.abund[j,1]],]
    out = specR[[specRinds3[1]]]$species.ordered })
  uniqDesignSpeciesNames = unique(designSpeciesNames);  nuniqDesignSpeciesNames = length(uniqDesignSpeciesNames)
  indRefSpeciesList = lapply(1:nuniqDesignSpeciesNames,function(j){
    vec=which(designSpeciesNames==uniqDesignSpeciesNames[j])
    out=list(species=uniqDesignSpeciesNames[j],index=vec)})
  indRefSpecies = sapply(1:nuniqDesignSpeciesNames,function(j){indRefSpeciesList[[j]]$index})
  nuniqRefSpecies = length(indRefSpecies)
  indOp = matrix(c(c(1,1+cumsum(num_ops_per_strains[1:(nmock-1)])),cumsum(num_ops_per_strains[1:nmock])),nrow=nmock,ncol=2)
} else if (substr(whichMock,start=1,stop=4)=="mSer"){ # The various Sereika et al. designs
  # Notes for D6322 strains:-
  #    Bacillus subtilis reclassified to spizizenii; ATCC 6633. specR
  #    Enterococcus faecalis strain B-537  - NCBI Reference Sequence: NZ_CP117970.1. sequence submitted by Zymo Research Corp. Feb 2023
  #    Escherichia coli
  designStrainNames = c("Bacillus subtilis B354 ATCC 6633 CP034943","Enterococcus faecalis B357 ATCC7080 GCF_028743535.1",
                        "Escherichia coli B1109", "Listeria monocytogenes B33116 ATCC19117",
                        "Pseudomonas aeruginosa B3509 ATCC15442","Salmonella enterica B4212 TA1536",
                        "Staphylococcus aureus B41012")
  nmock = length(designStrainNames)
  if (whichMock=="mSerF"){
    RA = c(13.2, 18.8,11.0,17.9,7.9,11.2,19.6)   # cell number.  Equal masses of genomic DNA. (== Opcounts/Number_of_operons_per_genome)
    RA = RA/sum(RA)
  } else {
    RA0 = c(13.2, 18.8,11.0,17.9,7.9,11.2,19.6)
    # Sub-sampling rates, ordered by designStrainNames order above
    #   Sub1:  0.1     1      1      1    1   0.5   0.2
    #   Sub2:  0.5    0.01   0.1     1    1    1     1
    #   Sub3:   1     0.02    1      1    1    1    0.02
    #   Sub4:   1      1     0.1     1    1    1    0.02
    
    
    if (whichSubMock==1){
      wts = c(0.1, 1, 1, 1, 1, 0.5, 0.2)
      RA1 = RA0*wts
      RA = RA1/sum(RA1)
    } else if (whichSubMock==2){
      wts = c(0.5, 0.01, 0.1, 1, 1, 1, 1)
      RA1 = RA0*wts
      RA = RA1/sum(RA1)
    }  else if (whichSubMock==3){
      wts = c( 1, 0.02, 1, 1, 1, 1, 0.02)
      RA1 = RA0*wts
      RA = RA1/sum(RA1)
    }  else if (whichSubMock==4){
      wts = c(1, 1, 0.1, 1, 1, 1, 0.02)
      RA1 = RA0*wts
      RA = RA1/sum(RA1)
    }
  }

  designSpeciesNames = sapply(1:length(designStrainNames),function(j){
    t1 = unlist(strsplit(designStrainNames[j],split=" "))[1:2]
    out=paste(t1,sep=" ",collapse=" ")})
  indOp = t(cbind(c(1,10),c(11,14),c(15,21),c(22,27),c(28,31),c(32,38),c(39,44)))
  num_ops_per_strains = indOp[,2] - indOp[,1] + 1
} else { # The various Srinivas et al. datasets.
  whichPair = as.numeric(substr(mstr,start=3,stop=3))
  if (s2c(mstr)[2] == "A"){
    designSpeciesNames = c("Bacillus subtilis", "Chromobacterium violaceum", "Enterococcus faecalis",
                           "Escherichia coli", "Halobacillus halophilus", "Haloferax volcanii", 
                           "Micrococcus luteus", "Pseudoalteromonas translucida","Pseudomonas fluorescens",
                           "Staphylococcus epidermidis")
    designStrainNames = designSpeciesNames
    genomeLengths = c(4214814,4750000,2739625,4639221,4170707,4018418,2501097,5021465,6511547,2497508)
    gM = rep(0.1,length(designSpeciesNames))
    nopsLratio = c(2.373,1.680,1.460,1.509,1.678,0.498,0.800,0.199,0.768,2.002)
    RAops = gM*nopsLratio;    RAops = RAops/sum(RAops)
    indOp = rbind(c(1,10),c(11,18),c(19,22),c(23,29),c(30,36),c(37,38),c(39,40),c(41,41),c(42,46),c(47,51))
    num_ops_per_strains = indOp[,2] - indOp[,1] + 1
    RAcells = gM/genomeLengths;  RAcells = RAcells/sum(RAcells)
    RA = RAcells
  } else if (s2c(mstr)[2] == "Z"){
    designSpeciesNames = c("Listeria monocytogenes", "Pseudomonas aeruginosa", "Bacillus subtilis",
                             "Escherichia coli", "Salmonella enterica", "Limosilactobacillus fermentum", 
                             "Enterococcus faecalis","Staphylococcus aureus")
    designStrainNames = designSpeciesNames
    genomeLengths = c(2880000,6770000,4076630,4773399,4780000,1796949,2900000,2872769)
    gM = c(0.898,0.0898,0.00898,0.000898,0.0000898,0.00000898,0.000000898,0.000000090)
    nopsLratio = c(2.083, 0.591, 2.453, 1.466, 1.464, 0.556, 1.379, 2.089)
    RAops = gM*nopsLratio;    RAops = RAops/sum(RAops)
    indOp = rbind(c(1,6),c(7,10),c(11,20),c(21,27),c(28,34),c(35,35),c(36,39),c(40,45)) 
    num_ops_per_strains = indOp[,2] - indOp[,1] + 1
    RAcells = gM/genomeLengths;  RAcells = RAcells/sum(RAcells)
    RA = RAcells
  } else {
    cat("No valid value for mstr. STOP. \n")
    stop()
  }
}
nmock = length(designStrainNames)

if (whichMock %in% c("mock50","mockKB")){
  iSpec = sapply(1:nmock,function(j){which(sapply(1:length(specR),function(k){specR[[k]]$species.ordered}) == designSpeciesNames[j])})
  
  if (whichMock=="mock50"){
    # Ensure compatibility of M (below) with that used to generate this mock microbiome.
    # Select a value of M to define the range of relative abundances from 1 to M - e.g. M=1000
    M = 1000
    RA = generate_RA(nmock,M)
  } else if (whichMock=="mockKB"){
    RA = RAkbg$Average
  } else {  # The whichMock == "mock10" case
    # From King etal's GFDB the relative abundances for the selected species are
    abund.Z = c(0.15, 0.0378, 0.0218, 0.0112, 0.00491,0.00259, 0.00185, 0.000705, 0.000046, 0.000042)
    RA = abund.Z/sum(abund.Z)
  }
  RA = as.vector(RA)
  outname = paste(dataset,"_",ERstring,"_characterisation.RData",sep="")
  if (substr(mstr,start=1, stop=3) %in% c("mSA","mSZ")){
    save(designStrainNames,specR,RA,iSpec,indOp,file=file.path(outRDatapath,outname))
  } else {
    save(designStrainNames,kbgood,ord.abKB,jord.abund,KB,specR,KBspecR,RA,iSpec,indOp,file=file.path(outRDatapath,outname))
  }
}    #    end   mock<nn> conditional    block

# If the kmer spectra for the ASVs and their associated reads have not been computed then call the 
# shell script  ../nanoporeSimulation/shellScripts/ASVkmerSpectra.sh  <whichSubunit>  <nthreads>
# noting that this may take an hour or so to complete if the microbiome dataset is large - e.g.
# greater than 100,000 reads. Generation of UMAP plots - see PART 4, ~ line 1350 - depends on this
# having been run for the particular amplicon (whichSubunit).
# Sample call is 
#    sbatch --mem=32GB --time=2:00:00 /vast/projects/rrn/nanoporeSimulation/shellScripts/ASVkmerSpectra.sh rrn 56
if(computeKmerSpectra){
  cat("\n Computing kmer spectra of all ASVs via a batch job dispatched at this stage of the code.\n")
  argstr1 = paste(" --mem=256GB --time=3:00:00 ",file.path(basepath,"ASVkmerSpectra.sh"), basepath,
                      whichMock,whichSubunit,whichCase,whichSubMock,whichPair,uppERstr,56,sep=" ")
  system2("sbatch",args=argstr1)
  cat("Code will be set to sleep for 15 minutes to allow computation of the kmer spectra.  Adjust as necessary. \n")
  Sys.sleep(900)
} else {
  cat("\n ASV kmer spectra assumed to have been previously computed. \n")
}

cat("Completed PART 0.\n")

# PART 1: Set up and call blastn on the dataset of interest (output fasta from denoising).
loThresh = ifelse(whichSubunit=="16S",1100,ifelse(whichSubunit=="23S",2100,3500))   # 2100 for 23S,   1400 for 16S
if (whichSubunit=="16S"){
  mmthresh = mmthresh16
  gapthresh = gapthresh16
} else if (whichSubunit=="ITS"){
  mmthresh = mmthreshITS
  gapthresh = gapthreshITS
} else if(whichSubunit=="23S"){
  mmthresh = mmthresh23
  gapthresh = gapthresh23
} else {
  mmthresh = mmthreshrrn
  gapthresh = gapthreshrrn
  
} 
denoise_indices_name = paste(stem1,"_",uppERstr,"_filtered_denoise_indices_out.txt", sep="")
denoise_names_name = paste(stem1,"_",uppERstr,"_filtered_denoise_names_out.txt", sep="")

DNpath = RADFastaPath
abind = 3
maxTargetSeqsA = 50
maxTargetSeqsB = 50
blastn_options_str1 = ' -outfmt "7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore score" '
blastn_options_str2A = paste('  -max_target_seqs ',maxTargetSeqsA,' -num_threads ',numcores,sep="",collapse="")
blastn_options_str2B = paste('  -max_target_seqs ',maxTargetSeqsB,' -num_threads ',numcores,sep="",collapse="")
blastn_call_db = paste(" -db ",file.path(blastDBpath,blastDBname),sep="")
blastn_call_query = paste(" -query ",file.path(RADFastaPath,ASVfastaName),sep="")
outnameA = paste("blastn_",stem1,"_",ERstring,"_",blastDBname.short,"_A.txt",sep="")
outnameB = paste("blastn_",stem1,"_",ERstring,"_",blastDBname.short,"_B.txt",sep="")
checkStr = paste("blastn_",stem1,"_",ERstring,"_",blastDBname.short,sep="")
btdir = dir(file.path(basepath,"text"),pattern="bl")
if (length(btdir)>0){
  ibt = which(sapply(1:length(btdir),function(j){t1 = str_locate(btdir[j], checkStr); out=!(is.na(t1[1,1]))}))
  if (length(ibt)>0){
    innameA = btdir[ibt[1]]   # This avoids the problem of knowing what date the text file was generated.
    innameB = btdir[ibt[2]]
  } else {
    innameA = paste("blastn_",stem1,"_",ERstring,"_",blastDBname.short,"_A.txt",sep="")
    innameB = paste("blastn_",stem1,"_",ERstring,"_",blastDBname.short,"_B.txt",sep="")
  }
} else {
  innameA = paste("blastn_",stem1,"_",ERstring,"_",blastDBname.short,"_A.txt",sep="")
  innameB = paste("blastn_",stem1,"_",ERstring,"_",blastDBname.short,"_B.txt",sep="")
}
blastn_call_outA = paste(" -out", file.path(basepath,"text",outnameA))
blastn_call_outB = paste(" -out", file.path(basepath,"text",outnameB))
argstrA = paste(blastn_call_db,blastn_call_query,blastn_options_str2A,blastn_options_str1,
                    blastn_call_outA,sep="")
argstrB = paste(blastn_call_db,blastn_call_query,blastn_options_str2B,blastn_call_outB,sep="")
# Following system2() call took 1020 secs. for 765 ASVs each of length about 4500 bases. With 8 threads
# it took 268 secs., with 16 threads took 141 secs..
if (doAlign){
  system2("/vast/projects/rrn/GROND/ncbi-blast-2.16.0+/bin/blastn",args=argstrA)
  system2("/vast/projects/rrn/GROND/ncbi-blast-2.16.0+/bin/blastn",args=argstrB)
}

cat("\n\n Completed Part 1. \n\n")

# PART 2: Outputs basic abundance calculation (genus, species, strain) plus ..
# Part 2.1: Parsing the A and B forms of output.  B form only provides the detailed name of the reference.
#           Key result is generation of Sobs which gives alignment details for each ASV.
#           This also includes filtering for unsatisfactory alignments (with ASVs potentially removed 
#           from further consideration)

if (doAlign){
  minAlen = ifelse(whichSubunit=="16S",1100,ifelse(whichSubunit=="ITS",250,ifelse(whichSubunit=="23S",2100,3800)))
  if (minAlen==3800 && substr(mstr,start=3,stop=3) %in% c("3","4")){minAlen=3300}
  A = readLines(con=file.path(basepath,"text",innameA))
  B = readLines(con=file.path(basepath,"text",innameB))
  kstartA = which(sapply(1:length(A),function(j){grep("Query:",A[j])})>0)
  kstartB = which(sapply(1:length(B),function(j){grep("Query=",B[j])})>0)
  kstartC = which(sapply(1:length(B),function(j){grep("^>",B[j])})>0)
  nASVs = length(kstartA)
  outname = paste("alignment_parsing_",dataset,"_",ERstring,".RData",sep="")
  save(minAlen,kstartA,kstartB,kstartC,A,B,file=file.path(basepath,"RData",outname))
} else {
  inname = paste("alignment_parsing_",dataset,"_",ERstring,".RData",sep="")
  load(file=file.path(basepath,"RData",inname))
  nASVs = length(kstartA)
}
# How many hits - derived from fileB. Note that number of hits is sometimes as small as 2
nhits = rep(0,nASVs)
for (k in 1:nASVs){
  ib = which.min(sapply(1:length(kstartC),function(j){d = kstartC[j]-kstartB[k]; out=ifelse(d<0,10000,d)}))
  delta = kstartC[ib]-kstartB[k]
  nhits[k] = delta-8
}
stats1 = data.frame(abund=rep(0,nASVs),ID=rep("",nASVs),accver=rep("",nASVs))
ab = rep(0,nASVs)
nhitsA = rep(0,nASVs)
S = vector("list",nASVs)
NS = vector("list",nASVs)
removedASVs = rep(0,10);  nlost = 0
kc = 0
for (k in 1:nASVs){
  diff = ifelse(k<nASVs,kstartA[k+1]-kstartA[k],length(A)-kstartA[nASVs]+1)
  ab[k] = as.numeric(unlist(strsplit(A[kstartA[k]],split="_"))[3])
  # ab[k] = as.numeric(unlist(strsplit(A[kstartA[k]],split="_"))[ifelse(whichDenoiser==1,2,4)])
  if (diff>(min(10,maxTargetSeqsA)+2)){
    nhitsA[k] = as.numeric(unlist(strsplit(A[kstartA[k]+3],split=" "))[2])
    numh = rep(0,nhitsA[k]);   chrh = rep("",nhitsA[k])
    stats2 = data.frame(ID=chrh,accver=chrh,PID=numh,Alen=numh,MM=numh,Gap=numh,Eval=numh,bitscore=numh,score=numh)
    for (k2 in 1:nhitsA[k]){
      line1 = unlist(strsplit(A[kstartA[k]+3+k2],"\t"))
      stats2[k2, 2:ncol(stats2)] = line1[c(2,3,4,5,6,11,12,13)]
      line2 = unlist(strsplit(B[kstartB[k]+5+k2],"  ")); i2 = which(nchar(line2)>1)
      stats2[k2,1] = line2[i2[1]]
    }
  } else { # Probably have zero hits
    nhitsA[k] = nhits[k]
    if (nhitsA[k]>0){
      numh = rep(0,nhitsA[k]);   chrh = rep("",nhitsA[k])
      stats2 = data.frame(ID=chrh,accver=chrh,PID=numh,Alen=numh,MM=numh,Gap=numh,Eval=numh,bitscore=numh,score=numh)
      for (k2 in 1:nhitsA[k]){
        line1 = unlist(strsplit(A[kstartA[k]+3+k2],"\t"))
        stats2[k2, 2:ncol(stats2)] = line1[c(2,3,4,5,6,11,12,13)]
        line2 = unlist(strsplit(B[kstartB[k]+5+k2],"  ")); i2 = which(nchar(line2)>1)
        stats2[k2,1] = line2[i2[1]]
      }
    } else {
      numh = 0;   chrh = "--"
      stats2 = data.frame(ID=chrh,accver=chrh,PID=numh,Alen=numh,MM=numh,Gap=numh,Eval=numh,bitscore=numh,score=numh)
    }
  }
  # Insert a filtering here that rejects unsatisfactorily aligned ASVs.
  # The only filter, initially (7Sept2023), will be on the basis of alignment length, stats2$Alen .
  # Note that bitscore could also be reasonably used.
  if (as.numeric(stats2[1,"Alen"]) > minAlen){
    kc = kc + 1
    S[[kc]] = stats2
  } else { # Record which of the ASVs are being excluded - with indexing based on the full set of ASVs
    nlost = nlost + 1
    removedASVs[nlost] = k
    NS[[nlost]] = stats2
  }
}
removedASVs = removedASVs[1:nlost]
cat("Number of ASVs removed because they have too limited alignment extent is ",nlost, " of ",nASVs,"\n")
iab = setdiff(1:length(ab),removedASVs)
ab = ab[iab]
nhits = nhits[iab]
nhitsA = nhitsA[iab]
if (whichSubunit == "16S"){
  S16 = S[1:kc]
  nASVs16 = nASVs  
} else if (whichSubunit == "ITS"){
  SITS = S[1:kc]
  nASVsITS = nASVs  
} else if (whichSubunit == "23S"){
  S23 = S[1:kc]
  nASVs23 = nASVs 
} else if (whichSubunit == "rrn"){
  Srrn = S[1:kc]
  nASVsrrn = nASVs 
} 
cat("Have removed ASVs ",removedASVs, " from further consideration (original ASV set indexing). \n")
nASVs = kc
Sobs = S[1:kc]
NS = NS[1:nlost]

decodeCrap = nchar(Sobs[[1]]$accver[1]) < 12
if (decodeCrap){
  # The database operons have been given simple, but not very informative, IDs - e.g. crap<nnn> -
  # due to difficulty in getting  makeblastdb  to recognise that there are no duplicates in the 
  # informative sequence names as stored in specR.  Need to load the file with specR in it from
  # the code that generated the fasta file used in creating the reference database. Then substitute
  # the IDs given in the list object S generated above by the corresponding informative IDs.
  inname = "specR_etal_grondDB.RData"   # This gives specR in which elements of operon name are separated by spaces.
  load(file=file.path("/vast/projects/rrn/GROND/output/RData",inname))
  OpNames = rep("",100000)
  opsCount = 0
  for (j in 1:length(specR)){
    nops = length(specR[[j]]$strainOps)
    OpNames[(opsCount+1):(opsCount+nops)] = specR[[j]]$strainOps
    opsCount = opsCount + nops
  }
  OpNames = OpNames[1:opsCount]
  for (ja in 1:nASVs){
    for (k in 1:nhits[ja]){
      id = Sobs[[ja]]$accver[k]
      if (nchar(id)<5){cat("ASV",ja,"  Hit",k," has id of length",nchar(id),"  ID is",id,"\n")}
      m = as.integer(substr(id,start=5,stop=nchar(id)))
      newID = OpNames[m]
      Sobs[[ja]]$ID[k] = newID
    }
  }
}

cat("Completion of Part 2.1 \n")

# Part 2.2: Operon, strain, species and genus relative abundance calculations and plot.

# First some initialising of generic (genus+species+strain+operon) material
genusNames = rep("",nASVs)
speciesNames = rep("",nASVs)
strainNames = rep("",nASVs)
strainNamesfull = rep("",nASVs)
operonNames = rep("",nASVs)
operonNamesfull = rep("",nASVs)
genusNamesPoor = rep("",nASVs)
speciesNamesPoor  = rep("",nASVs)
strainNamesPoor  = rep("",nASVs)
strainNamesfullPoor  = rep("",nASVs)
operonNamesPoor  = rep("",nASVs)
operonNamesfullPoor  = rep("",nASVs)


# Collect genus, species, strain and operon names of the ref DB to which each ASV best aligned.
#     Note that this code can be generalised to look at alignments that are almost as good 
#     as the best.
#        e.g. use whichRank =1 for best alignment, > 1 for others
whichRank = 1
for (k in 1:nASVs){
  # Check that the alignment meets quality criteria. If not assign to classification "other", and also
  # record the genus, species and strain details of that poor match.
  is.poor = (as.numeric(S[[k]]$MM[whichRank])>mmthresh) || (as.numeric(S[[k]]$Gap[whichRank])>gapthresh)
  if (is.poor){
    genusNames[k] = "other"
    speciesNames[k] = "otherg others"
    strainNames[k] = "othert"
    strainNamesfull[k] = "otherg others othert"
    operonNames[k] = "othert op1"
    operonNamesfull[k] = "otherg others othert op1"
    splitID = splitNamesSxy(Sobs[[k]]$ID[whichRank])
    genusNamesPoor[k] = splitID[1]
    speciesNamesPoor[k] = splitID[2]
    operonNamesPoor[k] = splitID[3]
    operonNamesfullPoor[k] = splitID[4]
  } else {
    splitID = splitNamesSxy(Sobs[[k]]$ID[whichRank])
    genusNames[k] = splitID[1]
    speciesNames[k] = splitID[2]
    strainNamesfull[k] = splitID[3]
    operonNamesfull[k] = splitID[4]
    specNamelength1 = nchar(unlist(strsplit(speciesNames[k],split=" "))[2])
    if (specNamelength1<4){
      specNamelengthR = specNamelength1
      jwR=1
      bitscore1 = as.numeric(Sobs[[k]]$bitscore[whichRank])
      while (specNamelengthR<4 && jwR<nhits[k]){
        # Check the alignment one rank lower for this ASV
        bitscore2 = as.numeric(Sobs[[k]]$bitscore[whichRank+1])
        if ((1-bitscore2/bitscore1)<0.001){
          splitID = splitNamesSxy(Sobs[[k]]$ID[jwR+1])
          genusNames[k] = splitID[1]
          speciesNames[k] = splitID[2]
          strainNamesfull[k] = splitID[3]
          operonNamesfull[k] = splitID[4]
          specNamelengthR = nchar(unlist(strsplit(speciesNames[k],split=" "))[2])
        }
        jwR = jwR+1
      }
      cat("For ASV ",k, " chose alignment ranked ",jwR,"\n")
    }
  }
} 
uniqOperonNamesfull = unique(operonNamesfull);    nuniqOperonNamesfull = length(uniqOperonNamesfull)
uniqStrainNamesfull = unique(strainNamesfull);    nuniqStrainNamesfull = length(uniqStrainNamesfull)
uniqSpeciesNames = unique(speciesNames);          nuniqSpeciesNames = length(uniqSpeciesNames)
uniqGenusNames = unique(genusNames);              nuniqGenusNames = length(uniqGenusNames)

cat("Observed number of unique genera, species, strains and operons ",nuniqGenusNames,nuniqSpeciesNames,
            nuniqStrainNamesfull,nuniqOperonNamesfull,"\n")


# Create specRfullC,
#fastaName = "rrnSpeciesDB_grond_refseq207full_C.fasta"
#specRfullC  = create_specRfullC(fastaName)
#nsp = length(specRfullC)
#save(specRfullC,file=file.path(basepath,"RData","specRfullC.RData"))
load(file.path(basepath,"RData","specRfullC.RData"))

cat("\nCompletion of Part 2.1: initialisation \n")

# Part 2.2: Compute abundances at different taxonomic resolutions and using a selected
# allocation criterion, and generate scatterplots of observed vs. design(expected) relative abundances.

# Case 1: Operon level         Allocation based on single top alignment.
# First, identify all operons occurring at top of list of ASVs and near-ASVs.
uniqOperonNamesfull = unique(operonNamesfull);     nuniqOperonNamesfull = length(uniqOperonNamesfull) 
zz = which(uniqOperonNamesfull=="")  # Need to deal with cases of zero counts for some ASVs, and hence an empty name 
if (length(zz)>0){uniqOperonNamesfull = uniqOperonNamesfull[setdiff(1:length(uniqOperonNamesfull),zz)]}
nuniqOperonfull = length(uniqOperonNamesfull)
operonCounts = rep(0,nuniqOperonfull) 
ASVtoOperonMap = rep(0,nASVs)
for (k in 1:nASVs){
  k3 = which(uniqOperonNamesfull==operonNamesfull[k])
  if (length(k3)>0){
    ASVtoOperonMap[k] = k3
    operonCounts[k3] = operonCounts[k3] + ab[k]
  }
}
OperonASVlist = lapply(1:nuniqOperonfull,function(k){which(ASVtoOperonMap==k)})

# Abundance Estimation at Operon Level.

#ops_per_species = vector("list",nuniqSpecies)
opsCounts_by_species  = vector("list",nuniqSpeciesNames)
opsNames_by_species  = vector("list",nuniqSpeciesNames)
uniqOperonCounts = sapply(1:nuniqOperonfull,function(jop){sum(ab[OperonASVlist[[jop]]])})
for (jst in 1:nuniqSpeciesNames){
  # Which OperonASVlist elements give ASV indices for operons of this strain?
  kndOp = which(sapply(1:nuniqOperonNamesfull,function(j){
    t1= unlist(strsplit(uniqOperonNamesfull[j],split=" "))
    species = paste(t1[1:2],sep=" ",collapse=" ")
    out = species == uniqSpeciesNames[jst]
  } ) )
  kndOp.ordered = kndOp[order(uniqOperonNamesfull[kndOp])]
  opsCounts_by_species[[jst]] = sapply(1:length(kndOp.ordered),function(j){uniqOperonCounts[kndOp.ordered[j]]})
  opsNames_by_species[[jst]] = sapply(1:length(kndOp.ordered),function(j){uniqOperonNamesfull[kndOp.ordered[j]]})
}

kord =  order(uniqOperonNamesfull)
X1 = cbind(uniqOperonNamesfull[kord],operonCounts[kord])
X2 = vector("list",nuniqStrainNamesfull)

Op.strains = sapply(1:nuniqOperonfull,function(j){OptoStrain(uniqOperonNamesfull[j])})
X3 = t(sapply(1:nuniqStrainNamesfull,function(j){
          i1 = which(Op.strains == uniqStrainNamesfull[j])
          thisCount = as.numeric(operonCounts[i1])
          out = c(uniqStrainNamesfull[j],paste(thisCount,sep="_",collapse="_"))
     })  )
iDtoO = sapply(1:nmock,function(j){t1=which(designStrainNames==uniqStrainNamesfull[j]); out=ifelse(length(t1)>0,t1,0)})
iOtoD = sapply(1:nuniqStrainNamesfull,function(j){t1=which(uniqStrainNamesfull==designStrainNames[j]); out=ifelse(length(t1)>0,t1,0)})


# Case  2: Strain level - by condensation of operon data.
uniqStrainNamesfull = unique(strainNamesfull)
zz = which(uniqStrainNamesfull=="")  # Need to deal with cases of zero counts for some ASVs, and hence an empty name 
if (length(zz)>0){uniqStrainNamesfull = uniqStrainNamesfull[setdiff(1:length(uniqStrainNamesfull),zz)]}
nuniqStrainfull = length(uniqStrainNamesfull)
strainCounts = rep(0,nuniqStrainfull) 
ASVtoStrainMap = rep(0,nASVs)
for (k in 1:nASVs){
  k3 = which(uniqStrainNamesfull==strainNamesfull[k])
  if (length(k3)>0){
    ASVtoStrainMap[k] = k3
    strainCounts[k3] = strainCounts[k3] + ab[k]
  }
}
StrainASVlist = lapply(1:nuniqStrainfull,function(k){which(ASVtoStrainMap==k)})
iobs = sapply(1:nmock,function(j){t1 = which(uniqStrainNamesfull == designStrainNames[j])
                                            out = ifelse(length(t1)==0,0,t1)})
# The index of iobs is for designStrainName, the value of iobs is the index in uniqStrainNamesfull of 
# the same name. If there was no observation of designStrainName[j] then iobs[j]==0
# We want obsProps to give zero for any designStrainNames that were not observed, and a 
# non-zero value for those that were observed.  It has no information on invalid reports of 
# an observation.
obsProps = rep(0,nmock)
for (j in 1:nmock){
  if (iobs[j]==0){
    obsProps[j] = 0
  } else {
    obsProps[j] = strainCounts[iobs[j]]/num_ops_per_strains[j]
  }
}
imisObs = setdiff(1:nuniqStrainNamesfull,iobs[which(iobs>0)])
# Construct Dstrain
if (length(imisObs)>0){
  # The following uses a species value for the number of operons of strains not in the designed
  # microbiome, or - if such is not available - the value 5.
  cat("Create Nopset giving the number of operons in each of these strains ('other' has 5) \n")
  Nopset = rep(0,length(imisObs))
  for (k in 1:length(imisObs)){
    specname = OptoSpecies(uniqStrainNamesfull[imisObs[k]])
    p = which(designSpeciesNames==specname)
    if (length(p)>0){
      m = indOp[p[1],2] - indOp[p[1],1] + 1
    } else {
      iz = which_specR(specname)
      if (length(iz)>0){
        m = length(specR[[iz]]$idb)
        if (m>17){
          m = 5
        }
      } else {m = 5}
    }
    Nopset[k] = m
  }
  cat("Invalid strains are \n")
  print(uniqStrainNamesfull[imisObs])
  make_DmisStrain = function(imisObs,Nopset){
    out = data.frame(obsStrain=uniqStrainNamesfull[imisObs],
                     obsCounts=strainCounts[imisObs],
                     obsNops= Nopset,
                     obsProps=strainCounts[imisObs]/Nopset)
  }
  
  nObsError = nuniqStrainNamesfull - length(iobs[which(iobs>0)])
  DmisStrain = make_DmisStrain(imisObs,4)
  Dstrain = data.frame(designStrain=c(designStrainNames,DmisStrain$obsStrain),
                                       obsProps=c(obsProps,DmisStrain$obsProps),
                                        designProps=c(RA,rep(0,length(imisObs))))
} else {
  Dstrain = data.frame(designStrain=designStrainNames,obsProps=obsProps,designProps=RA)
}
Dstrain$obsProps = Dstrain$obsProps /sum(Dstrain$obsProps )
print(Dstrain)
nr = length(Dstrain$designStrain)
nobs = nuniqStrainNamesfull

strainPlot = TRUE
plotName = paste("scatterplots_",dataset,"_",ERstring,".pdf",sep="")
pdf(file=file.path(basepath,"plots",plotName),paper="a4r")
nofitstrain = !strainPlot
if (strainPlot){
  ivalid = intersect(which(Dstrain$designProps>0.00001),which(Dstrain$obsProps>0.00001))
  if (length(ivalid)<2){
    cat("Fewer than 2 valid strains observed.  No plot generated for strains. \n")
    nofitstrain = TRUE
  } else {
    imisReported = which(Dstrain$designProps<0.0000001);  imissed = which(Dstrain$obsProps<0.0000001)
    if ((mstr=="SZ") || (mstr=="mockKB")){
      xv = log10(0.00001+Dstrain$designProps);     yv = log10(0.00001+Dstrain$obsProps);  nxv = length(xv)
      xva = xv[ivalid];   yva = yv[ivalid]   # Only used for coefficient calculation.
      plottitle = paste("Strain log10(Proportions)",dataset, ERstring, sep=" ")
      plot(xv, yv,main=plottitle, col="white",
           xlab="log10(Design proportions)",  ylab="log10(Observed Proportions)",
           sub=paste("(Correlation log(observed RA) vs. log(Design RA) ",round(100*cor(xva,yva))/100,")"),
           xlim= c(-5,3),   ylim = c(-5,0))
      for (j in 1:nr){points(xv[j],yv[j],pch=(j %% 25)+1,col=(j %% 24)+1)}
      legend(x="topright",legend=Dstrain$designStrain[ivalid], title="Valid",
             pch=(c(ivalid) %% 25)+1 , col=(c(ivalid) %% 24)+1,cex = 0.8,pt.cex=1.5)
      inotvalid = setdiff(1:nr,ivalid)
      if (length(inotvalid)>0){
        if (whichSubunit =="16S"){xnv = -1;  ynv = max(yv)} else {xnv=-4.9; ynv=0}
        legend(x=-1.1,y=-1.9,legend=Dstrain$designStrain[inotvalid],title="Missed or Invalid",
               pch=inotvalid %% 25+1 , col=inotvalid %% 24+1,cex = 0.8,pt.cex=1.5)
    #    legend(x=-4.8,y=0.2,legend=Dstrain$designStrain[inotvalid],title="Missed or Invalid",
    #           pch=inotvalid %% 25+1 , col=inotvalid %% 24+1,cex = 0.8,pt.cex=1.5)
      }
      fit.strain = lm(yva ~ xva)
      abline(reg=fit.strain,col="red")
    } else {
      xv = Dstrain$designProps;     yv = Dstrain$obsProps
      xva = xv[ivalid];   yva = yv[ivalid]   # Only used for coefficient calculation.
      plottitle = paste("Strain Proportions",dataset,ERstring,sep=" ")
      plot(xv, yv,main= plottitle, col="white",
           xlab="Design proportions",  ylab="Observed Proportions",
           sub=paste("(Correlation observed RA vs. Design RA ",round(100*cor(xva,yva))/100,")"),
           xlim= c(0,1.5*max(xv)),   ylim = c(0,max(yv)))
      for (j in 1:nr){points(xv[j],yv[j],pch=(j %% 25)+1,col=(j %% 24)+1)}
      legend(x=0.105,y=max(yv),legend=Dstrain$designStrain[ivalid],title="Valid",
             pch=(ivalid %% 25)+1 , col=(ivalid %% 24)+1)
      inotvalid = setdiff(1:nr,ivalid)
      if (length(inotvalid)>0){
        xnv = 0.01;  ynv = max(yv)
        legend(x=xnv,y=ynv,legend=Dstrain$designStrain[inotvalid],title="Missed or Invalid",
               pch=inotvalid %% 25+1 , col=inotvalid %% 24+1, cex = 1, pt.cex = 1.5) 
      }
      fit.strain = lm(yva ~ xva)
      abline(reg=fit.strain,col="red")
    }
  }    #   end  conditional (ivalid)   block
} else {fit.strain=NULL}

#######         Get species data directly from the strain data by condensation.       ########
AllSpecs = sapply(1:length(Dstrain$designStrain),function(j){
  t1=unlist(strsplit(Dstrain$designStrain[j],split=" "))[1:2]
  out = paste(t1,sep=" ",collapse=" ")})
uniqAllSpecs = unique(AllSpecs);             nuniqAllSpecs = length(uniqAllSpecs)
specObs = sapply(1:nuniqAllSpecs,function(j){
  iset = which(AllSpecs == uniqAllSpecs[j])
  out = sum(Dstrain$obsProps[iset])})
specObs.design = sapply(1:nuniqAllSpecs,function(j){
  iset = which(AllSpecs == uniqAllSpecs[j])
  out = sum(Dstrain$designProps[iset])})
Dspecies = data.frame(Species=uniqAllSpecs, obsProps=specObs, designProps=specObs.design)
SpeciesofStrains = sapply(1:nuniqStrainNamesfull,function(j){out=OptoSpecies(uniqStrainNamesfull[j])})
if (substr(mstr,start=1,stop=2) %in%c("SA","SZ")){SpeciesASVlist = StrainASVlist} else {
  StraintoSpecies = sapply(1:nuniqSpeciesNames,function(j){which(SpeciesofStrains==uniqSpeciesNames[j])})
  SpeciesASVlist = sapply(1:nuniqSpeciesNames,function(j){
                               vec=NULL
                               if (uniqSpeciesNames[j] == "otherg others"){out = NULL}else {
                                 nstrains = length(StraintoSpecies[[j]])
                                 cat(uniqSpeciesNames[j]," nstrains ",nstrains,"\n")
                                 for (k in 1:nstrains){
                                   vec = append(vec,StrainASVlist[[StraintoSpecies[[j]][k]]])
                                 }
                                 out = vec
                               }
  } )
}

speciesPlot = TRUE
if (speciesPlot){
  nr = length(Dspecies$Species)
  ivalid = intersect(which(Dspecies$designProps>0.00001),which(Dspecies$obsProps>0.00001))
  imisReported = which(Dspecies$designProps<0.0000001);  imissed = which(Dspecies$obsProps<0.0000001)
  if ((mstr=="SZ") || (mstr=="mockKB")){
    xv = log10(0.00001+Dspecies$designProps);     yv = log10(0.00001+Dspecies$obsProps);  nxv = length(xv)
    xva = xv[ivalid];   yva = yv[ivalid]   # Only used for coefficient calculation.
    plottitle = paste("Species log10(Proportions)",dataset, ERstring, sep=" ")
    plot(xv, yv,main=plottitle, col="white",
         xlab="log10(Design proportions)",  ylab="log10(Observed Proportions)",
         sub=paste("(Correlation log(observed RA) vs. log(Design RA) ",round(100*cor(xva,yva))/100,")"),
         xlim= c(-5,3),   ylim = c(-5,0))
    for (j in 1:nr){points(xv[j],yv[j],pch=(j %% 25)+1,col=(j %% 24)+1)}
    legend(x="topright",legend=Dspecies$Species[ivalid], title="Valid",
           pch=(c(ivalid) %% 25)+1 , col=(c(ivalid) %% 24)+1,cex = 0.9,pt.cex=1.5)
    inotvalid = setdiff(1:nr,ivalid)
    if (length(inotvalid)>0){
      if (whichSubunit =="16S"){xnv = -4.8;  ynv = 0} else {xnv=-4.9; ynv=0}
      legend(x=xnv,y=ynv,legend=Dspecies$Species[inotvalid],title="Missed or Invalid",
           pch=inotvalid %% 25+1 , col=inotvalid %% 24+1,cex = 0.8,pt.cex=1.5)
    }
    fit.species = lm(yva ~ xva)
    abline(reg=fit.species,col="red")
  } else {
    xv = Dspecies$designProps;     yv = Dspecies$obsProps
    xva = xv[ivalid];   yva = yv[ivalid]   # Only used for coefficient calculation.
    fit.species = lm(yva ~ xva)
    gradfit = summary(fit.species)$coefficients[2,1];  intcptfit = summary(fit.species)$coefficients[1,1]
    fitStr = paste("  LsqFit: Grad",round(100*gradfit)/100,"  Intcpt",round(100*intcptfit)/100,sep=" ")
    plottitle = paste("Species Proportions",dataset,ERstring,sep=" ")
    plot(xv, yv,main=plottitle, col="white",
         xlab="Design proportions",  ylab="Observed Proportions",
         sub=paste("(Correlation observed RA vs. Design RA ",round(100*cor(xva,yva))/100,fitStr,")",sep=" "),
         xlim= c(0,1.5*max(xv)),   ylim = c(0,max(yv)))
    for (j in 1:nr){points(xv[j],yv[j],pch=(j %% 25)+1,col=(j %% 24)+1)}
#    legend(x=0.105,y=max(yv),legend=Dspecies$Species[ivalid],title="Valid",
#           pch=(ivalid %% 25)+1 , col=(ivalid %% 24)+1)
    legend(x="bottomright",legend=Dspecies$Species[ivalid],title="Valid",
           pch=(ivalid %% 25)+1 , col=(ivalid %% 24)+1)
    inotvalid = setdiff(1:nr,ivalid)
    if (length(inotvalid)>0){
      xnv = 0.01;  ynv = max(yv)
      legend(x=xnv,y=ynv,legend=Dspecies$Species[inotvalid],title="Missed or Invalid",
           pch=inotvalid %% 25+1 , col=inotvalid %% 24+1, cex = 1, pt.cex = 1.5) 
    }
    fit.species = lm(yva ~ xva)
    abline(reg=fit.species,col="red")
  }
}

#######         Get genus data directly from the species data by condensation.       ########
AllGenus = sapply(1:length(Dspecies$Species),function(j){
  t1=unlist(strsplit(Dspecies$Species[j],split=" "))[1]
  out = t1})
uniqAllGenus = unique(AllGenus);             nuniqAllGenus = length(uniqAllGenus)
genusObs = sapply(1:nuniqAllGenus,function(j){
  iset = which(AllGenus == uniqAllGenus[j])
  out = sum(Dspecies$obsProps[iset])})
genusObs.design = sapply(1:nuniqAllGenus,function(j){
  iset = which(AllGenus == uniqAllGenus[j])
  out = sum(Dspecies$designProps[iset])})
Dgenus = data.frame(Genus=uniqAllGenus, obsProps=genusObs, designProps=genusObs.design)

genusPlot = TRUE
if (genusPlot){
  ivalid = intersect(which(Dgenus$designProps>0.00001),which(Dgenus$obsProps>0.00001))
  imisReported = which(Dgenus$Genus<0.0000001);  imissed = which(Dgenus$obsProps<0.0000001)
  if ((mstr=="SZ") || (mstr=="mockKB")){
    xv = log10(0.00001+Dgenus$designProps);     yv = log10(0.00001+Dgenus$obsProps);  nxv = length(xv)
    xva = xv[ivalid];   yva = yv[ivalid]   # Only used for coefficient calculation.
    plottitle = paste("Genus log10(Proportions)",dataset, ERstring, sep=" ")
    plot(xv, yv,main=plottitle, col="white",
         xlab="log10(Design proportions)",  ylab="log10(Observed Proportions)",
         sub=paste("(Correlation log(observed RA) vs. log(Design RA) ",round(100*cor(xva,yva))/100,")"),
         xlim= c(-5,3),   ylim = c(-5,0))
    for (j in 1:nr){points(xv[j],yv[j],pch=(j %% 25)+1,col=(j %% 24)+1)}
    legend(x="topright",legend=Dgenus$Genus[ivalid], title="Valid",
           pch=(c(ivalid) %% 25)+1 , col=(c(ivalid) %% 24)+1,cex = 0.9,pt.cex=1.5)
    inotvalid = setdiff(1:length(Dgenus$Genus),ivalid)
    if (length(inotvalid)>0){
      if (whichSubunit =="16S"){xnv = -0.5;  ynv = -1.5} else {xnv=-4.9; ynv=0}
      legend(x=-4.8, y=0,legend=Dgenus$Genus[inotvalid],title="Missed or Invalid",
           pch=inotvalid %% 25+1 , col=inotvalid %% 24+1,cex = 0.8,pt.cex=1.5)
    }
    fit.genus = lm(yva ~ xva)
    abline(reg=fit.genus,col="red")
  } else {
    xv = Dgenus$designProps;     yv = Dgenus$obsProps
    xva = xv[ivalid];   yva = yv[ivalid]   # Only used for coefficient calculation.
    plottitle = paste("Genus Proportions",dataset,ERstring,sep=" ")
    plot(xv, yv,main=plottitle, col="white",
         xlab="Design proportions",  ylab="Observed Proportions",
         sub=paste("(Correlation observed RA vs. Design RA ",round(100*cor(xva,yva))/100,")"),
         xlim= c(0,1.2*max(xv)),   ylim = c(0,max(yv)))
    for (j in 1:nr){points(xv[j],yv[j],pch=(j %% 25)+1,col=(j %% 24)+1)}
    legend(x="bottomright",legend=Dgenus$Genus[ivalid],title="Valid",
           pch=(ivalid %% 25)+1 , col=(ivalid %% 24)+1)
    inotvalid = setdiff(1:length(Dgenus$Genus),ivalid)
    if (length(inotvalid)>0){
      xnv = 0.01;  ynv = max(yv)
      legend(x=-4.8,y=0,legend=Dgenus$Genus[inotvalid],title="Missed or Invalid",
           pch=inotvalid %% 25+1 , col=inotvalid %% 24+1, cex = 1, pt.cex = 1.5) 
    }
    fit.genus = lm(yva ~ xva)
    abline(reg=fit.genus,col="red")
  }
}
 
dev.off()
 
outname=paste("blastn_species_Details_",dataset,"_",ERstring,".RData",sep="")
if (nofitstrain){
  save(uniqStrainNamesfull,StrainASVlist,strainCounts,
      uniqSpeciesNames, nofitstrain, fit.species,fit.genus,
       opsCounts_by_species,opsNames_by_species,uniqOperonCounts, 
        Dstrain, Dspecies, Dgenus, specR, specRfullC,Sobs,removedASVs,
     file=file.path(basepath,"RData",outname))
} else {
  save(uniqStrainNamesfull,StrainASVlist,strainCounts,
       uniqSpeciesNames, fit.strain, fit.species,fit.genus,
       opsCounts_by_species,opsNames_by_species,uniqOperonCounts, 
       Dstrain, Dspecies, Dgenus, specR, specRfullC,Sobs,removedASVs,
       file=file.path(basepath,"RData",outname))  
}

cat("\n Completed Part 2 - scatterplots \n\n")

# PART 3.   Generation of 2D UMAP Projections and Plots    

# Murrell group code, sequmap(), is available to compute 2D UMAP coordinates from the denoising output.
# However for large libraries I have found that this code takes too long to run (no parallelism, 100K 23S or rrn 
# reads do not complete within 48 hrs).  Hence an alternate R package, uwot, has been used which allows one to
# generate a UMAP model using a subset of the data and then embed the remaining data. This has been adopted here
# to reduce the computation for large libraries. 
# 
# Part 3.1: Import text files generated from Julia code call of denoise() followed by 
#           writing of ASV indices to text file.
#           Plot read error rate histogram for the fastq file processed by RAD.
#           If projJulia is TRUE Plot the UMAP data used in seqUMAP() of Murrell group, annotating strains.
#            - include plots that present only a single species, where that species 
#                has more than 1 strain present .
#            - also include a plot that has the less abundant strains not already 
#                presented in the single species plots.

stem2 = "filtered_denoise"
taskSet = c("indices","names","proj","templates")
currentrunID = dataset

# Jindices gives the indices into the filtered fastq file (that is fed to RAD) that
# are associated with each of the ASVs - so Jindices[37] is a string of multiple 
# indices separated by \t associated with ASV 37.  
# Jnames gives the headers of each of the fastq records in the filtered fastq file.
# P is the matrix of 2D UMAP coordinates of the ASVs and reads.
# Templ gives the sequences of each of the ASVs. 
if (projJulia){textfileSet = c(1,2,3,4)} else {textfileSet = c(1,2,4)} 
for (whichtask in textfileSet){
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
ASVreadIndices = sapply(1:nASV,function(j){as.numeric(unlist(strsplit(Jindices[j],split="\t")))})
kS = 6

lenThreshLow =   ifelse(whichSubunit=="16S",1100, ifelse(whichSubunit=="23S",2100,4300))
lenThreshHi =   ifelse(whichSubunit=="16S",1600, ifelse(whichSubunit=="23S",2950,5600))
if (prefiltplot){
  # Compute and plot read quality for pre-filtered files being processed by Julia code.
  inputName = paste(stem1,"_",ERstring,".fastq",sep="" ) 
  plotName = paste(stem1,"_",ERstring,".pdf",sep="" )
  #  inputName = paste(paste(whichMock,whichSubunit,CaseNumStr,"filtered",sep="_"),".fastq",sep="" ) # For SA<n> cases  
  if (nchar(whichMock)>3){ 
    # The mock microbiome dataset is generated with simulated reads and so the fastq sequence names of 
    # carry PID data. This is not the case for the mock microbiomes whose read data is real nanopore-
    # sequenced reads. The code of function 
    # 
     TT = plot_quality_histo(inputName, plotName,ThreshLow=lenThreshLow,ThreshHi=lenThreshHi)
    ifilt1 = TT[[1]];    T.df = TT[[2]]
  } else {
    cat("No prefilt plot is provided for real sequenced datasets.\n")
  }
}

# Compute and plot read quality for post-filtered files being processed by Julia code.

inputName = paste(stem1,"_",ERstring,"_filtered.fastq",sep="" )
plotName = paste(stem1,"_",ERstring,"_filtered.pdf",sep="" ) 
#  inputName = paste(paste(whichMock,whichSubunit,CaseNumStr,"filtered",sep="_"),".fastq",sep="" ) # For SA<n> cases  
iret = plot_quality_histo(inputName, plotName,ThreshLow=lenThreshLow,ThreshHi=lenThreshHi)


if (projJulia){ 
  # 2D UMAP projection directly from RAD(==Julia) run, but blastn run on this RAD run must be available.
  # Some datasets are strain level while others are species level.  
  # Strain level sets are "mock50", "mock59", "mockKB"
  # Species level sets are "mock10", "mockSA" "mockSZ", "mSerF"
  # Load from blastn_call_analysis_v02.R  
  #    uniqStrainNamesfull,StrainASVlist,strainCounts,
  #    opsCounts_by_strain,opsNames_by_strain,uniqOperonCounts
  plotname = paste("UMAP2D_Julia_",stem1,"_",ERstring,".pdf",sep="")
  pdf(file=file.path(basepath,"plots",plotname), paper="a4r")
  if (substr(mstr,start=1,stop=2) %in% c("SA","SZ","mS")){ # ,"mSerF","mSerS1","mSerS2","mSerS3","mSerS4")){
    # Load ab,uniqSpeciesNames,genusNames,SpeciesASVlist,speciesCounts,
    #    opsCounts_by_species,opsNames_by_species,uniqOperonCounts, 
    #      Dspecies, Dgenus, specR,specRfullC,Sobs,removedASVs
    inname=paste("blastn_species_Details_",dataset,"_",ERstring,".RData",sep="")
    load(file=file.path(basepath,"RData",inname))  
    nuniqSpeciesNames = length(uniqSpeciesNames)
    plottitle = paste("2D UMAP for ",dataset, ERstring," species labelled ASVs",sep=" ")
    plot.default(P[1,1:nASV],P[2,1:nASV],pch="+",col="white",
                 xlab="DIM 1",ylab="DIM 2", main=plottitle,
                 xlim = c(min(P[1,]),2*max(P[1,])))
    points(P[1,(nASV+1):ncol(P)],P[2,(nASV+1):ncol(P)],pch=19,col="yellow")
    legend(x="topright",legend=uniqSpeciesNames,pch=((1:nuniqSpeciesNames) %% 25) + 1,
           col=((1:nuniqSpeciesNames) %% 24) + 1)
    
    # Because some ASVs may have been removed as part of the blast analysis, and 
    # SpeciesASVlist has been calculated following such removal, there is a 
    # mis-match between the indexing in that list and that in the UMAP projections
    # file.  This needs to be corrected prior to labelling.
    v1 = 1:nASV       # where nASV is the number of ASVs returned by RAD
    if (removedASVs[1]>0){
      v1[removedASVs] = -1
      v2 = which(v1>0)
    } else v2 = v1
    for (js in 1:nuniqSpeciesNames){
      
      pset = SpeciesASVlist[[js]]
      #nrem = length(removedASVs)
      #if (removedASVs[1] > 0){  # instead of    if (nrem > 0)
      #  vec = SpeciesASVlist[[js]]
      #  jc = 0
      #  while (jc<nrem){
      #    iadj = which(vec>removedASVs[nrem-jc])
      #    vec[iadj] = vec[iadj] + 1
      #    jc = jc + 1
      #  } 
      #  pset = vec
      #}
      points(P[1,v2[pset]],P[2,v2[pset]],pch=(js %% 25) + 1,col=(js %% 24) + 1)
    } 
  }else { # We are dealing with the simulated reads datasets, which have a lot more strains.
    #  inname=paste("blastn_strainDetails_",dataset,".RData",sep="")
    inname=paste("blastn_species_Details_",dataset,"_",ERstring,".RData",sep="")
    load(file=file.path(basepath,"RData",inname))
    nuniqStrainNamesfull = length(uniqStrainNamesfull)
    plot.default(P[1,1:nASV],P[2,1:nASV],pch="+",col="white",
                 xlab="DIM 1",ylab="DIM 2", main=paste("2D UMAP for ",dataset, ERstring," strains labelling of ASVs",sep=" "),
                 xlim = c(min(P[1,]),2*max(P[1,])))
    points(P[1,(nASV+1):ncol(P)],P[2,(nASV+1):ncol(P)],pch=19,col="yellow")
    legend(x="topright",legend=uniqStrainNamesfull,pch=((1:nuniqStrainNamesfull) %% 25) + 1,col=((1:nuniqStrainNamesfull) %% 24) + 1,cex = 0.7)
    for (js in 1:nuniqStrainNamesfull){
      pset = StrainASVlist[[js]]
      points(P[1,pset],P[2,pset],pch=(js %% 25) + 1,col=(js %% 24) + 1)
    }
    # Now plot for selected subsets of the strains - just Bifidobacterium longum, or E.coli, of Bacteroides fragilis
    speciesNames = sapply(1:nuniqStrainNamesfull,function(j){
      t1 = unlist(strsplit(uniqStrainNamesfull[j],split=" "))[1:2]
      out=paste(t1,sep=" ",collapse=" ")})
    speciesSet = c("Bifidobacterium longum","Bacteroides fragilis","Escherichia coli",
                    "Bifidobacterium bifidum","OtherLowRA","Six Highest RA")
    strainSubsets = vector("list",10)
    strainSubsets[[1]] = which(speciesNames=="Bifidobacterium longum")
    strainSubsets[[2]] = which(speciesNames=="Bacteroides fragilis")
    strainSubsets[[3]] = which(speciesNames=="Escherichia coli")
    strainSubsets[[4]] = which(speciesNames=="Bifidobacterium bifidum")
    xx=NULL; for (j in 1:4){xx = append(xx,strainSubsets[[j]])}
    strainSubsets[[5]] = setdiff(7:nuniqStrainNamesfull,xx)
    strainSubsets[[6]] = 1:6
    for (k in 1:6){
      i1vec =  strainSubsets[[k]] 
      if (length(i1vec)>1){
        plot.default(P[1,1:nASV],P[2,1:nASV],pch="+",col="white",
                     xlab="DIM 1",ylab="DIM 2", main=paste("2D UMAP for ",speciesSet[k],ERstring," with identified strains labelling of ASVs",sep=" "),
                     xlim = c(min(P[1,]),2*max(P[1,])))
        for (js in 1:length(i1vec)){
          pset = StrainASVlist[[i1vec[js]]]
          for (jp in pset){
            points(P[1,nASV+ASVreadIndices[[jp]]],P[2,nASV+ASVreadIndices[[jp]]],pch=19,col="yellow")
          }
          points(P[1,pset],P[2,pset],pch=(i1vec[js] %% 25) + 1,col=(i1vec[js] %% 24) + 1)
        }
        legend(x="topright",legend=uniqStrainNamesfull[i1vec],
               pch=(i1vec %% 25) + 1,
               col=(i1vec %% 24) + 1)
      } else {
        cat("\n strainSubset",k," has fewer than 2 entries. \n")
      }
    }
  }
  dev.off()
}

cat("\n      Completed   Part 3.1    \n")


# Part 3.2: Generation of 2D UMAP plots using uwot and potentially, sub-sampling and embedding.
#           Sub-sampling the read sets of ASVs with a large number of associated reads.
#           UMAP 2D plot with ASVs, training reads, and embedded reads identified
#
cat("\n Part 3.2 commenced. \n ")   
#plotname = paste("uwot-derived_UMAP2D_subsample_embed_",dataset,"_",ERstring,".pdf",sep="")
#pdf(file=file.path(basepath,"plots",plotname),paper="a4")
#cat("\n plotname : ", plotname,"\n")

# Load pre-computed kmer-based edit distances between ASVs + reads for the relevant dataset.
#CAVEAT: kmSpecNormed is an object holding kmer spectra, not kmer-based edit distances.
inname1 = paste("kmSpecNormed_",whichMock,"_",whichSubunit,"_Case_",whichCase,"_",ERstring,".RData",sep="")
load(file=file.path(basepath,"RData",inname1))

# Get key information to allow identification of ASVs with species.
# We are dealing with the simulated reads datasets, which have a lot more strains.
#  inname=paste("blastn_strainDetails_",dataset,".RData",sep="")
inname=paste("blastn_species_Details_",dataset,"_",ERstring,".RData",sep="")
load(file=file.path(basepath,"RData",inname))
nuniqStrainNamesfull = length(uniqStrainNamesfull)

# Multistrain species visualisation.

if (whichMock %in% c("mock50","mockKB")) {
  plotname = paste("uwot-derived_UMAP2D_rawManhattanMethod_selectedSpeciesmultiStrain_",dataset,"_",ERstring,".pdf",sep="")
  pdf(file=file.path(basepath, "plots",plotname),paper="a4r")
  cat("Full dataset UMAP2D Computation to be undertaken, using umap2. \n")
  cat("\n Starting plotting based on uwot calls for full read set. \n plotname ",plotname,"\n")
  
  
  # Now plot for selected subsets of the strains - just Bifidobacterium longum, or E.coli, or Bacteroides fragilis
  speciesNames = sapply(1:nuniqStrainNamesfull,function(j){
    t1 = unlist(strsplit(uniqStrainNamesfull[j],split=" "))[1:2]
    out=paste(t1,sep=" ",collapse=" ")})
  
  shortStrainNames = sapply(1:nuniqStrainNamesfull,function(j){
    t1 = unlist(strsplit(uniqStrainNamesfull[j],split=" "))
    out=t1[length(t1)]
  })
  
  strainSubsetLabel = c("Bifidobacterium longum","Bacteroides fragilis","Escherichia coli",
                        "Bifidobacterium bifidum","Other LowRA","Top 6 RA")
  strainSubsets = vector("list",10)
  strainSubsets[[1]] = which(speciesNames=="Bifidobacterium longum")
  strainSubsets[[2]] = which(speciesNames=="Bacteroides fragilis")
  strainSubsets[[3]] = which(speciesNames=="Escherichia coli")
  strainSubsets[[4]] = which(speciesNames=="Bifidobacterium bifidum")
  xx=NULL; for (j in 1:4){xx = append(xx,strainSubsets[[j]])}
  strainSubsets[[5]] = setdiff(7:nuniqStrainNamesfull,xx)                   # Other lower RA species
  strainSubsets[[6]] = 1:6                                                  # The 6 highest abundance strains
  Subsetproj.train = vector("list",6)
  DistMatsAll = vector("list",4)
  for (k in 1:4){
    strainSubsets[[k]];   nStSub = length(strainSubsets[[k]])   
    # strainSubsets[[k]] indexes the strains in the mock microbiome whose species name is, for instance (k==1), "Bifidobacterium longum" 
    ikmStrainSet = NULL
    ikmStrainSubsets = vector("list",nStSub)
    iStrainsASVSet = NULL;  countiSA = 0
    ik1 = strainSubsets[[k]]   #  ik1 is a set of indices into speciesNames, and also into uniqStrainNamesfull.
    for (k1 in 1:nStSub){
      cat("Processing strain ",uniqStrainNamesfull[strainSubsets[[k]]][k1],"\n")
      Ik1 = NULL
      ik2 = StrainASVlist[[ik1[k1]]]
      countiSA = countiSA + length(ik2)   # Counting the ASVs for this subset of strains
      iStrainsASVSet = append(iStrainsASVSet,ik2)   # Getting the indices of the ASVs for this strain subset
      for (k2 in 1:length(ik2)){
        ikmStrainSet = append(ikmStrainSet, ASVreadIndices[[ik2[k2]]])
        Ik1 = append(Ik1,ASVreadIndices[[ik2[k2]]])
      }
      ikmStrainSubsets[[k1]] = Ik1
    }
    iStrainsASVSet = iStrainsASVSet[1:countiSA];   niSA = length(iStrainsASVSet)
    # There are generally far more reads for each ASV than we can afford to compute 
    # distances for.  Hence we will sample a small number - e.g. 10 - per ASV.
    # So for Bifidobacterium longum there are 6 strains and 22 ASVs.  So if we have 10
    # reads sampled for each ASV there will be 220 reads.  Thus the distance measures
    # would be being computed for 22*(10+1)=242 rrn-length sequences.  More work needed
    # to get the code for this sorted out!
    # Ideally we would order the sequences as (ASV1, reads-of-ASV1, ASV2,reads-of-ASV2,...)
    sampsize = 20
    iseqs = NULL
    iASVs = NULL
    iASVlist = vector("list",nStSub)
    iASVreadslist = vector("list",nStSub)
    ireadsSubsampASVs = NULL
    jacount = 1
    for (k1 in 1:nStSub){      
      ik2 = StrainASVlist[[ik1[k1]]]
      iASVlist[[k1]] = ik2
      iASVvec = NULL
      t1 = NULL 
      for (k2 in 1:length(ik2)){
        iseqs = append(iseqs,ik2[k2])  # Add the ASV index
        iASVs = append(iASVs,jacount);  iASVvec = append(iASVvec,jacount)
        numASVReads = min(sampsize,length(ASVreadIndices[[ik2[k2]]]))
        isamp = sample(1:length(ASVreadIndices[[ik2[k2]]]),numASVReads)
        t1 = append(t1,jacount+1:numASVReads)
        jacount = jacount + 1 + numASVReads
        iseqs = append(iseqs, nASV+ASVreadIndices[[ik2[k2]]][isamp])
      }
      iASVlist[[k1]] = iASVvec 
      iASVreadslist[[k1]] = t1
    }
    # Now compute the set of distances between the operons of this species.
    
  #  editDistabs = function(sp1, sp2, k){
      # Computes the approximate Levenshtein distance between sequences
      # having length-nmormed kmer spectra sp1, sp2.
      #    sp1, sp2     the kmer spectrum for kmers of length  k  of sequences 1 and 2
      #    Lseq1, Lseq2 the length in bases of sequences 1 and 2.
      #  29 June 2025                                              [cjw]
 #     tot = sum(sapply(1:length(sp1),function(j){out=abs(sp1[j] - sp2[j])}))
 #     out = tot/(2*k)
 #   }
    
    EDkmerPar = function(i,allSpec,kS){
      # allSpec is a matrix of 4^kS rows, each column being a length-normed kmer spectrum
      # of a read for k-mer length kS.
      # Returns a vector of ncol(allSpecs) edit distances.
      # 29 June 2025                                               [cjw]
      sp1 = allSpec[,i]
      out = sapply((i+1):ncol(allSpec),function(j){sp2=allSpec[,j]; out=editDistabs(sp1,sp2,kS)})
    }
    
    
    
    iOpsSpec = kmSpecNormed[,iseqs]
    iOpsSpec.df = as.data.frame(t(iOpsSpec))
    
    n3 = ncol(iOpsSpec)  
    j3count = 0
    iOSdists = NULL   # OS from Observations x Samples where Observations are kmer spectra elements
    distveclist = mclapply(1:(n3-1),EDkmerPar,iOpsSpec,kS,mc.cores=numcores)
    distvec = NULL
    for (j3 in 1:length(distveclist)){
      distvec = append(distvec,unlist(distveclist[[j3]]))
    }
    iOSdists = distvec
      
    include2 = FALSE
    if (include2){
      for (j3 in 1:(n3-1)){
        #distvec = mclapply((j3+1):n3,EDkmerPar,iOpsSpec,kS,mc.cores=numcores) 
        distvec = NULL
        sp1 = iOpsSpec[,j3]
        for (i in (j3+1):n3){
          sp2 = iOpsSpec[,i]
          dd = editDistabs(sp1,sp2,kS)
          distvec = append(distvec,dd)
        }
        iOSdists = append(iOSdists,unlist(distvec))
        #    cat(j3,length(distvec),"\n")
        j3count = j3count + n3-j3
        #     rm(distvec)
      }     
    }     #   end    include2    conditional    block
    # Form an upper triangular matrix from iOSdists
    iOSdistmat = matrix(0, nrow = n3, ncol = n3)
    iOSdistmatu = matrix(0, nrow = n3, ncol = n3);  iOSdistmatl = iOSdistmatu
    iOSdistmatu[upper.tri(iOSdistmatu, diag = FALSE)] = iOSdists
    iOSdistmatl = t(iOSdistmatu)
    # iOSdistmatl[lower.tri(iOSdistmatl, diag = FALSE)] = iOSdists
    iOSdistmat = iOSdistmatl + iOSdistmatu
    image(t(iOSdistmat)/max(iOSdistmat))
    DistMatsAll[[k]] =  round(5000*iOSdistmat) 
 
    
    cat("About to call umap2 at line 1576 \n")
    # Am using a dataframe as input which consists of observations on 4^kS variables, the variables being 
    # kS-mers - 4^kS of them.  This is analogous to the demonstration material using the iris dataset.
    Subsetproj.train[[k]] = umap(round(5000*iOpsSpec.df), ret_model = FALSE, 
                                 metric="manhattan",n_neighbors=5, nn_method="nndescent") 
    Subsetproj.train[[k]] = umap(as.dist(iOSdistmat), ret_model = FALSE, 
                                 n_neighbors=5, nn_method="nndescent") 
for (k in 1:4){
    shortSN = shortStrainNames[strainSubsets[[k]]] 
    if (length(shortSN)<7){
      shortSNstr = paste(shortSN,sep=" ",collapse=" ")
    } else {
      shortSNstr = paste(paste(shortSN[1:5],sep=" ",collapse=" "),"...",sep="")
    }
    plot(Subsetproj.train[[k]][iASVs,1],Subsetproj.train[[k]][iASVs,2],xlab="DIM1",ylab="DIM2",
           main = paste("Species",strainSubsetLabel[k],"  Strains",shortSNstr,sep=" "),
            xlim = c(min(Subsetproj.train[[k]][,1]), 1.8*max(Subsetproj.train[[k]][,1])),
            pch=".",col="black")
    for (ja in 1:length(iASVlist)){
      points(Subsetproj.train[[k]][iASVlist[[ja]],1],Subsetproj.train[[k]][iASVlist[[ja]],2],pch=ja,col=ja+1, cex=1.8)
      points(Subsetproj.train[[k]][iASVreadslist[[ja]],1],Subsetproj.train[[k]][iASVreadslist[[ja]],2],pch=ja,col=ja+1, cex=0.5)
    }
    legend(x="bottomright",legend=shortSN,pch=1:length(iASVlist),col=1:length(iASVlist)+1)
    cat("Completed call of umap at line 1547.\n")
   
  }    #    end     k     loop
  dev.off()
  outname = "multistrain_distMat_umapCoordSets.RData"
  save(Subsetproj.train,DistMatsAll,strainSubsetLabel,file=file.path(basepath,"RData",outname))
  # load(file=file.path(basepath,"RData",outname))
}      #    end    conditional on being simulated data ("mock")    block



quit(save="no")


threshASVreads = 30
iumASV = vector("list",nASV)   # The set of indices of reads from each ASV that are to be 
                               # used in UMAP computation
iumASV = lapply(1:nASV,function(j){
  len = length(ASVreadIndices[[j]])
  if (len > threshASVreads){
    isamp  = sample(1:len,threshASVreads) 
    out = isamp 
  } else {
    out = 1:len
  }
})
#kSN = kmSpecNormed

# There is a need to deal with removed ASVs.  The following code is effective in the proJulia-conditioned
# plotting further on in this code.  Something equivalent is needed in what follows more immediately here.
#     
##    Because some ASVs may have been removed as part of the blast analysis, and 
##    SpeciesASVlist has been calculated following such removal, there is a 
##    mis-match between the indexing in that list and that in the UMAP projections
##    file.  This needs to be corrected prior to labelling.
##       
##       pset = SpeciesASVlist[[js]]
##       nrem = length(removedASVs)
##       if (removedASVs[1]> 0){  # instead of    if (nrem > 0)
##         vec = SpeciesASVlist[[js]]
##         jc = 0
##         while (jc<nrem){
##           iadj = which(vec>removedASVs[nrem-jc])
##           vec[iadj] = vec[iadj] + 1
##           jc = jc + 1
##         } 
##         pset = vec
##       }
##       points(P[1,pset],P[2,pset],pch=(js %% 25) + 1,col=(js %% 24) + 1)

colcount1=nASV;   colcount2 = nASV;  colcount3 = 0     # kmSpecNormed includes the ASV spectra in columns 1:nASV
inU = matrix(0,nrow=4^kS,ncol=(threshASVreads+1)*nASV)
inU[,1:nASV] = kmSpecNormed[,1:nASV]
inUset = matrix(0,nrow=nASV,ncol=2)
inotU = matrix(0,nrow=4^kS,ncol=length(Jnames))
innotUset = matrix(0,nrow=length(Jnames),ncol=2)
for (j in 1:nASV){
  nc = length(iumASV[[j]])
  #  inU[,(colcount1+1):(colcount1+nc)] = kmSpecNormed[,colcount2+iumASV[[j]]]
  inU[,(colcount1+1):(colcount1+nc)] = kmSpecNormed[,ASVreadIndices[[j]][iumASV[[j]]]]
  inUset[j,] = c(colcount1+1,colcount1+nc) 
  notc = length(ASVreadIndices[[j]]) - nc
  if (notc>0){
    ivec = setdiff(1:length(ASVreadIndices[[j]]),iumASV[[j]])
    inotU[,(colcount3+1):(colcount3+notc)] = kmSpecNormed[,nASV+setdiff(ASVreadIndices[[j]],ASVreadIndices[[j]][iumASV[[j]]])]
    innotUset[j,] = c(colcount3+1,colcount3+notc)
  }
  colcount1 = colcount1+nc
  colcount2 = colcount2 + length(ASVreadIndices[[j]])
  colcount3 = colcount3+notc
}
inU = inU[,1:colcount1]
inotU = inotU[,1:colcount3];    innotUset = innotUset[1:colcount3,] 

#proj = umap(inU, config = umap.defaults,method = "naive",n_neighbors=5,n_components=2,
#            metric="manhattan")

# Now construct the matrix whose columns are the kmer spectra of the reads indexed by iumASV.
# This will be the input to the call of umap.
# Also construct the complementary matrix whose columns are those not included in the umap input.
# This will be the input into umap.predict.
# Key function for edit distance calculations is
##     EDkmerPar = function(i,allSpec,kS){
##       allSpec is a matrix of 4^kS rows, each column being a length-normed kmer spectrum
##       of a read for k-mer length of kS.
##       Returns a vector of ncol(allSpecs) edit distances.
cat("Commencing computation of pairwise edit distances via kmer-approximation.\n")
cat("   First for all ASVs + no more than 30 reads per ASV. \n")
tic()
n1 = ncol(inU);  jcount = 0
inUdists = rep(0,n1*(n1-1) %/% 2)
for (j1 in 1:(n1-1)){
  distvec = mclapply((j1+1):n1,EDkmerPar,inU,kS,mc.cores=numcores) 
  inUdists = append(inUdists,distvec)
  jcount = jcount + n1-j1
}
toc()

cat("   Second for all reads not already considered. \n")
tic()
n2 = ncol(inotU);  j2count = 0
inotUdists = rep(0,n2*(n2-1) %/% 2)
for (j2 in 1:(n2-1)){
  distvec = mclapply((j2+1):n2,EDkmerPar,inotU,kS,mc.cores=numcores) 
  inotUdists = append(inotUdists,distvec)
  j2count = j2count + n2-j2
}
toc()

cat("\n Ready to start calls of umap from the uwot package. \n")
# Using uwot package
cat("Train umap \n")
tic()
proj.train = umap2(inUdists, ret_model = TRUE, nn_method="nndescent")
toc()
cat("Embed into umap \n")
tic()
proj.embed = umap_transform(t(colSums(inotU)),proj.train)
toc()
totkmerspec = ncol(kmSpecNormed)
P = t(proj.train$embedding)
Pem = t(proj.embed)
AllASVreadIndices = NULL
for (j in 1:nASV){AllASVreadIndices = append(AllASVreadIndices,unlist(ASVreadIndices[[j]])) }
readsnotASV = setdiff(1:(ncol(kmSpecNormed)-nASV+1),AllASVreadIndices)
proj.notASVreads = umap_transform(t(colSums(kmSpecNormed[,readsnotASV])),proj.train)
Pnot = t(proj.notASVreads)
plottitle = paste("UMAP (embedding) ",stem1, ERstring,sep=" ")


v1 = 1:nASV       # where nASV is the number of ASVs returned by RAD, nASVs is the number of retained ASVs.
if (removedASVs[1]>0){
  v1[removedASVs] = -1
  v2 = which(v1>0)
} else v2 = v1


xv = unlist(as.vector(P[1,v2]));     yv = unlist(as.vector(P[2,v2]))

par(mai=c(1.82,1.42,0.82,0.42))  # default is c(1.02,0.82,0.82,0.42)
par(omi=c(0.5,0.5,0,0))
magAxis = 1.8;  magLabel = 1.8;  magSub = 1.2
mySubstr = paste("Red: train reads (",ncol(P),");  Black: ASVs (",length(Jindices),");",
                 "   Green embedded reads (",ncol(Pem),");  Purple not ASV associated (",length(readsnotASV),")",sep=" ")
plot(xv,yv,pch="o",col="white",main=plottitle,xlab="Dim1",ylab="Dim2",
     cex.lab=magLabel, cex.axis=magAxis, cex.sub=magSub,
     xlim=c(min(xv),max(xv) + (max(xv)-min(xv))/2))
mtext(mySubstr, side=1, line=5, adj=0.5, cex=1.5, col="black",outer=FALSE)

points(Pem[1,],Pem[2,],pch="o",col="green",cex=0.4)                                                  # Plot UMAP2D embedded points
points(as.vector(P[1,(nASV+1):ncol(P)]),as.vector(P[2,(nASV+1):ncol(P)]),pch=19,col="red",cex=0.5)   # Plot UMAP2D training read points
points(as.vector(P[1,v2]),as.vector(P[2,v2]),pch="+",col="black",cex=0.9)                    # Plot UMAP2D ASV points
points(as.vector(Pnot[1,]),as.vector(Pnot[2,]),pch=18,col="#FF00FF")

for (js in 1:nuniqStrainNamesfull){
  # Plot, with distinctive markers for each strain, the UMAP2D ASVs for each strain
  # inU[,inUset[pset[1],1]:inUset[pset[1],2]]
  pset = StrainASVlist[[js]]
  points(P[1,v2[pset]],P[2,v2[pset]],pch=(js %% 25) + 1,col=(js %% 24) + 1,cex=1.5)
}
legend(x="topright",legend=uniqStrainNamesfull, cex=(1-0.005*nuniqStrainNamesfull),
       pch=((1:nuniqStrainNamesfull) %% 25) + 1,col=((1:nuniqStrainNamesfull) %% 24) + 1)

outname5 = paste(stem1,ERstring,"umap_projections_embedding.RData",sep="_")
save(proj.train,proj.embed,proj.notASVreads,file=file.path(basepath,"RData",outname5))

dev.off()

plotname = paste("uwot-derived_UMAP2D_allReads_",dataset,"_",ERstring,".pdf",sep="")
pdf(file=file.path(basepath, "plots",plotname),paper="a4r")
if (fullReadset){  
  cat("Full dataset UMAP2D Computation to be undertaken. \n")
  cat("\n Starting plotting based on uwot calls for full read set. plotname ",plotname,"\n")
  tic()
  proj.all = umap2(t(colSums(kmSpecNormed)), ret_model = TRUE, nn_method="nndescent")
  toc()
  cat("Completed uwot call on full kmerSpecNormed.\n")
  P = t(proj.all$embedding)
  plottitle = paste("UMAP all",stem1,ERstring,sep=" ")
  xv = unlist(as.vector(P[1,1:nASV]));     yv = unlist(as.vector(P[2,1:nASV]))   
  plot(xv,yv,pch="o",col="white",main=plottitle,xlab="Dim1",ylab="Dim2",
       sub=paste("Red: reads (",ncol(P),");  Black: ASVs (",length(Jindices),")",sep=" "),
       xlim=c(min(xv),max(xv) + 3*(max(xv)-min(xv))/5))
  points(as.vector(P[1,(nASV+1):ncol(P)]),as.vector(P[2,(nASV+1):ncol(P)]),pch=19,col="red",cex=0.5)
  points(as.vector(P[1,1:nASV]),as.vector(P[2,1:nASV]),pch="+",col="black",cex=0.8)
  legend(x="topright",legend=uniqSpeciesNames,pch=((1:nuniqSpeciesNames) %% 25) + 1,
         col=((1:nuniqSpeciesNames) %% 24) + 1)
  
  # Because some ASVs may have been removed as part of the blast analysis, and 
  # SpeciesASVlist has been calculated following such removal, there is a 
  # mis-match between the indexing in that list and that in the UMAP projections
  # file.  This needs to be corrected prior to labelling.
  v1 = 1:nASV       # where nASV is the number of ASVs returned by RAD
  if (removedASVs[1]>0){
    v1[removedASVs] = -1
    v2 = which(v1>0)
  } else v2 = v1
  for (js in 1:nuniqSpeciesNames){
    
    pset = SpeciesASVlist[[js]]
    #nrem = length(removedASVs)
    #if (removedASVs[1] > 0){  # instead of    if (nrem > 0)
    #  vec = SpeciesASVlist[[js]]
    #  jc = 0
    #  while (jc<nrem){
    #    iadj = which(vec>removedASVs[nrem-jc])
    #    vec[iadj] = vec[iadj] + 1
    #    jc = jc + 1
    #  } 
    #  pset = vec
    #}
    points(P[1,v2[pset]],P[2,v2[pset]],pch=(js %% 25) + 1,col=(js %% 24) + 1, cex=1.2)
  } 
  dev.off()
}

cat("\n\n      COMPLETED   PART  3    \n\n")

quit(save="no")


# PART 4: Gather data on ASV "purity", a measure of how close an ASV is to being associated with
#           a single operon. First form ASVPurityTables that, for each ASV, identify the origin of 
#           each read associated with this ASV - the specific operon - and how many reads are from 
#           each such operon identified. This leads to tables of operon counts for each strain 
#           represented by the reads associated with that ASV.
#           IMPORTANT:  ASV purity is only able to be determined from simulated reads datasets, which
#                       have sequence headers that contain operon-level identification.
if (whichMock %in% c("mockKB","mock50")){
  ASVreadIndices = sapply(1:nASV,function(j){as.numeric(unlist(strsplit(Jindices[j],split="\t")))})
  ASVtoOp = sapply(1:nASV,function(j){out = iret[ASVreadIndices[[j]]]})     
  # ASVtoOp description:
  #     For all ASVs gives the index into the input fastq file for each read associated with a single ASV
  # Is there a single strain or several strains associated with a particular ASV?
  ASVreadStrainNames = sapply(1:nASV,function(j){
    nops = length(T.df[ASVtoOp[[j]],"operon"])
    strNames = rep("",nops)
    for (jop in 1:nops){
      t1 = unlist(strsplit(T.df[ASVtoOp[[j]][jop],"operon"],split="[_]"))
      strNames[jop] = paste(t1[1:(length(t1)-1)],sep = "_",collapse = "_")
    }
    out = strNames  })
  ASVOpIndices = sapply(1:nASV,function(j){out=T.df[ASVtoOp[[j]],"operon"]})
  # Do any of the ASVs have more than 1 strain of read associated with them? If so which ones?
  multiStrainASVs = rep(0,nASV)
  jcount = 0
  for (ja in 1:nASV){
    if (length(unique(ASVreadStrainNames[[ja]]))>1){ 
      jcount = jcount + 1
      multiStrainASVs[jcount] = ja
      cat("\nASV ",ja, " has multiple strains associated with it.\n")
      print(unique(ASVreadStrainNames[[ja]]))
    }
  }
  multiStrainASVs = multiStrainASVs[1:jcount]
  nmultiSA = jcount
  cat("Number of ASVs having reads from multiple strains associated with them is",nmultiSA, "\n")
  
  # First deal with purity for those ASVs having a single strain.
  oneStrainSet = setdiff(1:nASV,multiStrainASVs);  noneSS = length(oneStrainSet)
  ASVPurityTables = vector("list",nASV)
  ASVPurityTables[oneStrainSet] = lapply(oneStrainSet,function(js){
    out=list(strainTables=ASVOpCount(ASVOpIndices[[js]]),
             strainNames=unique(ASVreadStrainNames[[js]]))
  })
  # Now deal with multi strain ASVs.
  for (js in multiStrainASVs){
    allStrains = unique(ASVreadStrainNames[[js]]);  nallStrains = length(allStrains)
    singleStrainTables = vector("list",nallStrains)
    for (k in 1:nallStrains){
      i1vec = which(ASVreadStrainNames[[js]]== allStrains[k])
      singleStrainTables[[k]] = ASVOpCount(ASVOpIndices[[js]][i1vec])
    }
    ASVPurityTables[[js]] = list(strainTables=singleStrainTables, strainNames=allStrains)
  }
  
  Purity =   lapply(1:nASV, function(j){
    if (j %in% multiStrainASVs){
      # Sum the opCounts for each strain and select that strain with the largest opCount.
      # Random choice if ties on opCount.
      nallStrains = length(ASVPurityTables[[j]]$strainTables)
      strainOpCounts = sapply(1:nallStrains,function(k){sum(ASVPurityTables[[j]]$strainTables[[k]]$opCount)})
      istr = which.max(strainOpCounts)
      bestStrain=istr
      bestCount=strainOpCounts[istr]
      totCounts = sum(strainOpCounts)
      strainPurity = bestCount/totCounts
      ordering = order(strainOpCounts,decreasing=TRUE)
      bestStrain = ASVPurityTables[[j]]$strainNames[istr]
      otherStrains = setdiff(ASVPurityTables[[j]]$strainNames,bestStrain)
    } else {
      strainOpCounts = ASVPurityTables[[j]]$strainTables$opCount
      istr = which.max(strainOpCounts)
      bestCount=sum(strainOpCounts)
      totCounts = bestCount
      strainPurity = bestCount/totCounts
      bestStrain = ASVPurityTables[[j]]$strainNames[1]
      otherStrains = NULL
      ordering = 0
    }
    out = list(strainPurity=strainPurity,dominantStrain=bestStrain, otherStrains=otherStrains,
               totCount=totCounts,ordering=ordering)  
  })
  
  sP = sapply(1:nASV,function(j){Purity[[j]]$strainPurity})
  dS = sapply(1:nASV,function(j){Purity[[j]]$dominantStrain})
  tC = sapply(1:nASV,function(j){Purity[[j]]$totCount})
  ASVpurity = data.frame(strainPurity=sP, dominantStrain=dS, totCounts=tC)
  
  cat("Total count of reads associated with ASVs is ",sum(ASVpurity$totCount), " of", nfiltr," filtered reads. \n")
  cat(" Count of reads associated with ASVs that are not 100% pure is ", sum(ASVpurity$totCount[which(ASVpurity$strainPurity<1)])," \n")
  impure = which(ASVpurity$strainPurity<1)
  print(cbind(impure,ASVpurity[impure,]))
  
  # Now identify all strains occurring in impure ASVs and the read counts allocated to each.  
  # Return as a table with columns headed    Name   Actual_Count   Count_as_Dominant
  
  actualCount=rep(0,nASV);   strainName=rep("",nASV);  asDomCount=rep(0,length(impure));  DomName = rep("",length(impure))
  jnames = 0;  jDomCount = 0
  for (ja in 1:length(impure)){
    theseNames = ASVPurityTables[[impure[ja]]]$strainNames;       ntNames = length(theseNames)
    for (j1 in 1:ntNames){
      i1 = which(strainName==theseNames[j1])
      if (length(i1)==0){ # This is a new strainName
        strainName[jnames+1] = theseNames[j1]
        actualCount[jnames+1] = actualCount[jnames+1] + sum(ASVPurityTables[[impure[ja]]]$strainTables[[j1]]$opCount)
        jnames = jnames+1
      } else {  # There should not be any name duplication, so length(i1)==1 
        actualCount[i1] = actualCount[i1] + sum(ASVPurityTables[[impure[ja]]]$strainTables[[j1]]$opCount)
        strainName[i1] = theseNames[j1]
      }
    }
    i2 = which(DomName ==Purity[[impure[ja]]]$dominantStrain)  # Also == theseNames[1]
    if (length(i2)==0){# This is a new dominant strain name
      DomName[jDomCount+1] = Purity[[impure[ja]]]$dominantStrain
      asDomCount[jDomCount+1] = asDomCount[jDomCount+1] + Purity[[impure[ja]]]$totCount
      jDomCount = jDomCount + 1
    } else {
      asDomCount[i2] = asDomCount[i2] + Purity[[impure[ja]]]$totCount
      #    DomName[i2] = Purity[[ja]]$dominantStrain
    }
  }
  
  impureCounts = list(allStrains=data.frame(Name=strainName[1:jnames],  Actual_Count=actualCount[1:jnames]), 
                      dominantStrains=data.frame(DominantNames=DomName[1:jDomCount], Count_as_Dominant=asDomCount[1:jDomCount]))
  
  for (ja in 1:length(impure)){
    cat("ASV ",impure[ja],"Purity Table \n")
    print(ASVPurityTables[[impure[ja]]])
  }
  
  notDominant = setdiff(impureCounts[[1]]$Name,impureCounts[[2]]$DominantNames)
  print(notDominant)
  nrow(impureCounts[[2]])
  cat("There are ", nrow(impureCounts[[1]])," strains that appear in impure ASVs. \n")
  cat("There are ",nrow(impureCounts[[2]]), "strains that appear in impure ASVs as dominant strains for at least one ASV.\n")
  cat("There are ",length(notDominant), "strains that do not appear in any ASV as dominant but do occur in at least one impure ASV. \n")
  
  # Which, if any, of the strains that are not dominant in any impure ASV do not have any (necessarily pure) ASV with
  # which they are associated?  The read count for such a strain will be lowered by its association with an impure
  # ASV.
  # Any strain that, if it occurs in any impureASV, always appears as the dominant strain will have its read count
  # incorrectly high.
  
  # Identify strains from the notDominant class that are associated with 1 or more pure ASVs.
  
  # Divert the output from the following loop to file. It is the key textual data on
  # the analysis of non-dominant strains in impure ASVs.
  outname = paste("nonDominant_reads_by_ASVs_",dataset,".txt",sep="")
  sink(file.path(basepath,"text",outname))
  nnD = length(notDominant)
  strainName = rep("",nnD)
  strainTotCount = rep(0,nnD)
  cat("Names and Counts for non-Dominant strains from impure ASVs dataset ",dataset,"\n")
  inotDom = sapply(1:nnD,function(j){which(impureCounts[[1]]$Name == notDominant[j])})
  print(impureCounts[[1]][inotDom,])
  cat("\n Names and Counts for Dominant strains from impure ASVs dataset ",dataset,"\n")
  iDom = setdiff(1:nrow(impureCounts[[1]]),inotDom)
  print(impureCounts[[1]][iDom,])
  missedStrains = matrix(0,nrow=nnD,ncol=2);  colnames(missedStrains) = c("designIndex", "Counts")
  # Need to convert all notDominant names into space-separated elements for matching.
  for (jnd in 1:nnD){
    cat("\n\n\nConsidering strain ",notDominant[jnd],"\n")
    # Is there at least one dominant strain that is identical to this strain?
    idom = which(strainNamesfull == gsub("_"," ",notDominant[jnd]))
    if (length(idom) == 0){
      # This strain has no ASV giving best alignment to this strain. So the 
      # strain is either missed or invalid.
      idesign = which(designStrainNames == gsub("_"," ",notDominant[jnd]))
      if (length(idesign)==0){
        # This strain is invalid. Leave its counts belong to the impure ASV to which 
        # this strain had been associated.
        cat("Strain ", notDominant[jnd], " is not a valid strain. \n")
        cat("  - its count remains with the dominant strain of the ASV with which this strain was associated.\n")
      } else {
        # This strain has been missed in the standard analysis.
        missedStrains[jnd,1] = idesign;   missedStrains[jnd,2] = impureCounts[[1]][inotDom[jnd],"Actual_Count"]
        cat("Strain ",notDominant[jnd], "is in the design, but was missed by the standard analysis.\n")
        cat("A pseudoASV can be created and the counts for that ASV are ",missedStrains[jnd,2],"\n" )
      }
    } else {
      # This strain is identical to one which is dominant in one or more ASVs - it is observed.
      # Add its count to strainNamesfull[idom].
      strainCounts[idom] = strainCounts[idom] + impureCounts[[1]][inotDom[jnd],"Actual_Count"]
      cat("Strain ",notDominant[jnd],"is identical to one which is dominant in one or more ASVs - it is observed. \n")
      cat("  - its count is added to the count of ",impureCounts[[1]][inotDom[jnd],"Actual_Count"]," from standard processing.\n")
    }
  }    #     end   jnd    loop
  sink()
  
  # For each observed operon identify the ASVs having that operons as their best alignment.
  # Provided by OperonASVlist[[jop]] for each operon, indexed by jop in 1:nuniqOperonNamesfull.
  # Then, for each ASV, jASV, in that list length(ASVreadIndices[[OperonASVlist[[jop]][jASV]]]) is 
  # the count of reads associated with that ASV.
  # Hence count for the operon is the sum of such read counts.
  
  OpCounts = rep(0,nuniqOperonNamesfull)
  
  for (jop in 1:nuniqOperonNamesfull){
    aset = OperonASVlist[[jop]]
    OpCounts[jop] = sum(sapply(aset,function(ja){length(ASVreadIndices[[ja]])}))
  }
  
  # Form a list of read counts for each observed operon of each strain.
  countStr = rep("",nuniqStrainNamesfull)
  OpStrains = sapply(1:length(operonNamesfull),function(j){OptoStrain(operonNamesfull[j])})
  OpsASVs = sapply(1:nuniqStrainNamesfull,function(j){which(OpStrains == uniqStrainNamesfull[j])})
  for (js in 1:nuniqStrainNamesfull){
    nops = length(unique(operonNamesfull[poolOps[[js]]]))
    countOps = rep(0,nops)
    for (jop in 1:nops){
      kv = which(operonNamesfull[OpsASVs[[js]]] == unique(operonNamesfull[OpsASVs[[js]]])[jop])
      countOps[jop] = sum(sapply(kv,function(kvj){length(ASVreadIndices[[unlist(OpsASVs[[js]])[kvj]]])}))
    }
    countStr[js] = paste(countOps[1:nops],sep="_",collapse="_")
  }
  strainOpsCounts = data.frame(strain=uniqStrainNamesfull, readCounts=countStr)
  print(strainOpsCounts)
  
  # Now estimate cellular strain counts using mean operon counts for each strain.
  strainCount = rep(0,nuniqStrainNamesfull)
  for (jstr in 1:nuniqStrainNamesfull){
    upp = indOp[jstr,2];   low = indOp[jstr,1]
    strainCount[jstr] = sum(OpCounts[low:upp])/(upp-low+1)
  }
  
  outname4 = paste(dataset,"_ASVtoOps_3November.RData",sep="")
  # load(file.path(basepath,outname3))
  save(ASVreadIndices,ASVPurity,ASVPurityTables,ASVOpIndices,ASVOpCount,Counts,OpCounts,strainCount,
       file=file.path(basepath,dataset,"RData",outname4))
} else {
  cat("Part 4 is only applicable to simulated reads mock microbiomes. Hence it has not been processed for dataset ",dataset,"\n")
}
cat("\n\n      Completed   Part 4    \n\n")





#############################################################################################
#############################################################################################
#########################    REMOVED CODE 22 November     ###################################
#############################################################################################
#############################################################################################

get_mode = function(x) {
  uniq_x = unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}



genus_in_specRfullC = function(genus){
  which(sapply(1:length(specRfullC),function(j){
    t1=unlist(strsplit(specRfullC[[j]]$species.ordered,split="_"))[1]}) == genus)}

genus_in_specR = function(genus){
  which(sapply(1:length(specR),function(j){
    t1=unlist(strsplit(specR[[j]]$species.ordered,split=" "))[1]}) == genus)}


editDistabs = function(sp1, sp2, k){
  # Computes the approximate Levenshtein distance between sequences
  # having length-nmormed kmer spectra sp1, sp2.
  #    sp1, sp2     the kmer spectrum for kmers of length  k  of sequences 1 and 2
  #    Lseq1, Lseq2 the length in bases of sequences 1 and 2.
  #  29 June 2025                                              [cjw]
  tot = sum(sapply(1:length(sp1),function(j){out=abs(sp1[j] - sp2[j])}))
  out = tot/(2*k)
}

EDkmerPar = function(i,allSpec,kS){
  # allSpec is a matrix of 4^kS rows, each column being a length-normed kmer spectrum
  # of a read for k-mer length of kS.
  # Returns a vector of ncol(allSpecs) edit distances.
  # 29 June 2025                                               [cjw]
  sp1 = allSpec[,i]
  out = sapply(1:ncol(allSpec),function(j){editDistabs(sp1,allSpec[,j],kS)})
}

readQscore = function(i,frags.Q,start_trim){
  t1 = as.character(frags.Q[[i]])
  t1c = unlist(strsplit(t1,split=""))
  t1.Q = sapply(1:nchar(t1),function(j){as.integer(charToRaw(t1c[j]))-33})
  t1.Qmean = -10*log10(mean(sapply(start_trim:nchar(t1),function(j){10^(-t1.Q[j]/10)})))
  t1.Qstd = -10*log10(sd(sapply(start_trim:nchar(t1),function(j){10^(-t1.Q[j]/10)})))
  t1.meanQ = mean(t1.Q)
  out = t1.Qmean       #  c(t1.Qmean,t1.Qstd,t1.meanQ)
}

BR_Qmeans = function(i,fastqPath,start_trim=61, makeplot2=TRUE){
  opnum=i 
  fastqName = paste("reads_",whichSubunit,"_Op_",opnum,".fastq",sep="")
  stem = unlist(strsplit(fastqName,split="[.]"))[1]
  inname1 = paste(stem,"_metadata_Pidfilt_Op_",opnum,".RData",sep="")
  load(file=file.path(basepath,"outputBR/RData",inname1))  # This loads  ifilt1  for this opnum
  pat = paste("Op_",opnum,".fastq",sep="",collapse="")
  F0 = readFastq(dirPath = fastqPath, pattern=pat)
  F = F0[ifilt1]
  # Compute read Qmean for all F. Written for the case of no parallelism, though
  # the function was written for parallel (mclapply) use.
  cat("opnum:  ",opnum,"  length(ifilt1): ",length(ifilt1),"\n")
  Qmean = rep(0,length(ifilt1))
  for (jf in 1:length(ifilt1)){
    Qmean[jf] = readQscore(jf,quality(F),start_trim)
  }
  fit = lm(y ~ x, data = data.frame(x=T.df[ifilt1,"PID"],y=Qmean))
  grad = round(1000*summary(fit)$coefficients[2,1])/1000  
  gradProb = round(1000*summary(fit)$coefficients[2,4])/1000
  fitStr = paste("(Linear fit: gradient ",grad," Prob(>|t|) ",gradProb,")",sep=" ")
  if (makeplot2){
    plotname2 = paste(stem,"_smoothScatter_Op_",opnum,".pdf",sep="")
    pdf(file=file.path(basepath,"outputBR/plots",plotname2),paper="a4")
    par(mfrow=c(1,1))
    plot(T.df[ifilt1,"PID"],Qmean,xlab="PID",ylab="Qmean", main="mock_23S_Op1")
    palette <- hcl.colors(30, palette = "inferno")
    smoothScatter(x=T.df[ifilt1,"PID"],y=Qmean,colramp = colorRampPalette(palette),
                  xlab="PID", main=paste("mock,  23S,  Operon ", opnum, sep=""), sub=fitStr)
    dev.off()
  }
  outname3 = paste(stem,"_Qmean.RData",sep="")
  save(Qmean,file=file.path(basepath,"outputBR/RData",outname3))
  out = Qmean
}


BRfastq_filtPID = function(i,fastqPath,whichSubunit,makeplot1=TRUE){ 
  opnum=i  
  stem = paste("reads_",whichSubunit,"_Op_",opnum,sep="")           #  unlist(strsplit(fastqName,split="[.]"))[1]
  readHeadersName = paste("reads_",whichSubunit,"_Op_",opnum,"_headers.txt",sep="")
  inpath = file.path(fastqPath,readHeadersName)
  X = readLines(con=file.path(basepath,"outputBR/text",readHeadersName))
  igood = which(sapply(1:length(X), function(j){nchar(X[j])< 180}))  
  newname = paste(stem,"_Qstripped.txt",sep="")
  writeLines(X[igood],con=file.path(basepath,"outputBR/text",newname))
  T = read.delim(file=file.path(basepath,"outputBR/text",newname),sep=" ")
  nr = nrow(T)
  colnames(T) = c("ID.ont","ID.mock","ReadLength","FragLength","ReadPID")
  print(T[1,])
  
  # Now parse ID.mock, ReadLength, FragLength, ReadPID to extract the key data.
  splitChar = c("",",","=","=","=")
  cname = c("ID.ont","ID.mock","ReadLength","FragLength","ReadPID")
  IDvec = rep("",nr);   Strand=rep(-9,nr)  
  ReadLen = rep(0,nr);  FragLen = rep(0,nr);   PID = rep(0.0,nr)
  for (jr in 1:nr){
    for (jc in 2:5){
      t1 = unlist(strsplit(T[jr,cname[jc]],split=splitChar[jc]))
      if (length(t1)>1){
        if (jc==2){
          cstrand = substr(t1[2],start=1,stop=1)
          Strand[jr] = ifelse(cstrand=="+",1,-1)
        } else {
          t11 = t1[2]
          if (jc==3){
            ReadLen[jr] = as.integer(t11)
          } else if (jc==4){
            FragLen[jr] = as.integer(t11)
          } else {
            PID[jr] = as.numeric(as.numeric(substr(t11,start=1,stop=nchar(t11)-1)))
          }
        }
      } else {
        cat("Improper split character for operon",opnum," row",jr, " column",jc,"\n")
        cat("   T[jr,cname[jc]]:  ",T[jr,cname[jc]]," split character:  ",splitChar[jc],"\n")
      }
    }
  }
  T.df = data.frame(ID=T[1:nr,1], strand=Strand, readLength=ReadLen, fragLength=FragLen,PID=PID)
  ifilt1 = which(sapply(1:length(T.df$PID),function(j){T.df$PID[j]>98.75}))  
  cat("Number of retained reads after initial filtering (PID=98.75) is ",length(ifilt1),"\n")
  if (makeplot1){
    plotname1 = paste(stem,"_histos_PID_Qmean_Op_",opnum,".pdf",sep="")
    pdf(file=file.path(basepath,"outputBR/plots",plotname1),paper="a4")
    par(mfrow=c(2,2))
    hist(T.df[,"PID"],xlab="PID",main=paste("Badread-simulated mock Op",opnum," PID",sep=" ") )
    hist(log10(100.001-T.df[,"PID"]),xlab="log10(100.001 - PID)",
         main=paste("Badread-simulated mock Op",opnum," PID",sep=" ")) #   which(T.df$PID<98.75)
    breaks1 = 100 - seq(98.74,100,by = 0.02)
    hist(1-T.df[ifilt1,"PID"]/100,breaks=breaks1/100,
         xlab="Fraction Identical",main=paste("Badread-simulated mock Op",opnum," PID",sep=" "))
    hist(log10(1.00001-T.df[ifilt1,"PID"]/100),
         xlab="log10(1.00001 - PID/100) (= log10(mismatch rate))",main=paste("Badread-simulated mock Op",opnum," PID",sep=" "))
    dev.off()
  }
  outname1 = paste(stem,"_metadata_Pidfilt_Op_",opnum,".RData",sep="")
  save(T.df,ifilt1,file=file.path(basepath,"outputBR/RData",outname1))
  out = list(Tdf = T.df, PIDfilt=ifilt1)
}


#######################################################################
# For Bifidobacterium bifidum strains we have the following:- 
#   For the 8 ASVs of the 3 strains
ivec = c(1,22,43,64,85,97,107,128)
print(round(DistMatsAll[[4]][ivec,ivec]))
# Splitting into the individual strains:-
ivecStr = vector("list",3)
ivecStr[[1]] = c(1,22);  ivecStr[[2]] = c(43,64,85,97);   ivecStr[[3]] = c(107,128)
round(5000*iOSdistmat[c(3,4,5),c(3,4,5,6)]/(2*kS-1))
for (j in 1:3){
  cat("Strain PRL2010");  print(round(DistMatsAll[[4]][ivecStr[[j]],ivecStr[[j]]]))
}

# Considering the reads now.  Take ASV372 of S17, which has 9 associated reads.
ivecR = 97:106
print(round(DistMatsAll[[4]][ivecR,ivecR]))


#######################################################################
#######################################################################
# Checking, using iris data, use of different forms of input to (uwot) umap.
nI = nrow(iris)
ns = 50
isamp = sample(1:nI,ns)
irisDist = matrix(1e-10,ns,ns)
for (j1 in 1:(ns-1)){
  for (j2 in (j1+1):ns){
    irisDist[j1,j2] = sum(sapply(1:4,function(j){abs(iris[isamp[j1],j]-iris[isamp[j2],j])}))
    irisDist[j2,j1] = irisDist[j1,j2]
  }
}
nn = 3
Umetric1 = "manhattan";  Umetric2 = "sumAbsDiff"
Iumap1 = umap(iris[isamp,1:4], ret_model = FALSE, 
              metric="manhattan",n_neighbors=nn) 
Iumap2 = umap(as.dist(irisDist), ret_model = FALSE, n_neighbors=nn)

par(mfrow=c(1,2))
plot(Iumap1[,1],Iumap1[,2],xlab="DIM1",ylab="DIM2",
     main = paste("Iris Data - Raw",nn,Umetric1,sep=" "),pch="o",col="black")
plot(Iumap2[,1],Iumap2[,2],xlab="DIM1",ylab="DIM2",
     main = paste("Iris Data - DistMat",nn,Umetric2,sep=" "),pch="+",col="blue")



for (ja in 1:length(iASVlist)){
  points(Subsetproj.train[[k]][iASVlist[[ja]],1],Subsetproj.train[[k]][iASVlist[[ja]],2],pch=ja,col=ja+1, cex=1.8)
  points(Subsetproj.train[[k]][iASVreadslist[[ja]],1],Subsetproj.train[[k]][iASVreadslist[[ja]],2],pch=ja,col=ja+1, cex=0.5)
}
legend(x="bottomright",legend=shortSN,pch=1:length(iASVlist),col=1:length(iASVlist)+1)


#  Subsetproj.train[[k]] = umap(round(5000*iOpsSpec.df), ret_model = FALSE, 
#                               metric="manhattan",n_neighbors=5, nn_method="nndescent") 
#Subsetproj.train[[k]] = umap(as.dist(iOSdistmat), ret_model = FALSE, 
#                             n_neighbors=5, nn_method="nndescent") 



