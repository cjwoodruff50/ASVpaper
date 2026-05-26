# This Rscript creates the design of the microbiome, extracts from the reference DB fasta 
# files corresponding to the operon or amplicon (16S and 23S rRNA genes or the rrn operon) sequences 
# of the bacteria in the mock microbiome, sets up for a call of the badread code that 
# simulates ONT nanopore-sequenced reads, and finally submits a set of batch jobs that 
# generate simulated ONT nanopore sequencer reads of each operon of the designed mock microbiome.
# The amplicon library for mock microbiomes based on the GROND refseq207full database and
# the King et al. (2019) proposed standard human healthy gut microbiome (referred to here as KBGF
# (Knowledge Base Gut Feeling) or, more commonly, as KB.
#
# source("/vast/projects/rrn/RscriptsArchive/current/make_mock_operon_store.R") 
#
# 15 December 2025                                                                    [cjw]
args = commandArgs(trailingOnly=TRUE)
#######################################################################################
##########################                               ##############################
##########################         INITIALISATIONS       ##############################
##########################                               ##############################
#######################################################################################
slurm = TRUE
if (length(args>0)){
  basepath = args[1]
  whichMock = args[2]
  whichSubunit = args[3]
  whichCase = args[4]
  identMean = args[5]
  identSD = args[6]
  maxopnum = args[7]
  nmock = args[8]
  maxTotfrags = as.numeric(args[9])
  numcores = as.numeric(args[10])
} else {
  basepath = "/vast/projects/rrn/ASVtest"
  whichMock = "mockKB"
  whichSubunit = "23S"
  whichCase = "11"
  identMean = 26
  identSD = 4
  maxopnum = 291  # 291 for mockKB (59 strains), 261 for mock50 (50 strains)
  nmock = 59
  maxTotfrags = 50000  
  numcores = 4
}

cat("Basepath  is ",basepath,"\n")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
packages <- c("ShortRead", "parallel","tictoc","readxl")    # Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))


#######################################################################################
##########################                               ##############################
##########################           FUNCTIONS           ##############################
##########################                               ##############################
#######################################################################################
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

linKBtospecR = function(inname,constructing=FALSE){
  # Creates key objects and mappings for 
  #    1. Identifying which KB strains have entries in specR;
  #    2. Ordering KB strains by relative abundance;
  #    3. Mapping abundance-ordered KB strains to their corresponding specR entries.
  #    4. Determining strain operons for any KB strain in specR.
  # 22 November 2025 
  # Load crapOpNames, opDB.species, opDB.strains, OpNames, ops_per_strain_each_spec, specR
  inname = "specR_etal_grondDB.RData"
  load(file.path(basepath,"RData",inname))
  
  # Get KB and KB abundance data to guide rank-ordering on KB abundance data of the strains being used 
  # in mock microbiomes.
  gfkb.name = "KGFDB.xlsx"  #  "King CH etal SuppMat S4 Table GutFeeling.KB PLoSOne.0206484.s008.xlsx"
  KB = read_excel(file.path(basepath,"GROND",gfkb.name))
  abundKB.name = "KGFDB_Abundance_Tables.xlsx"
  AbundKB = read_excel(file.path(basepath,"GROND",abundKB.name))
  
  if (constructing){
    # Build the mapping from KB to specR.  Doing this is a manual job that is assisted by the code in this 
    # conditional block.If the GROND strain database is changed this construction process may need updating.
    KBspecR = matrix(0,nrow=164,ncol=3)
    colnames(KBspecR) = c("ind.specR","opLow","opHi")
    # Manually repeat the following for j =1:164, saving the work after each value of j.
    outnameK = "KBspecR_matrix.RData"
    load(file.path(basepath,"RData",outnameK))
    # j=22
    # print(KB[j,c(3,8)])
    # iz = which_specR("Escherichia coli")
    # print(specR[[iz]]$strainOps)
    # KBspecR[j,] = c(0,-1,-1)
    # KBspecR[j,] = c(iz,1,2)
    # save(KBspecR,file=file.path(basepath,"RData",outnameK))
  }     #     end     conditional   constructing    block
  
  innameK = "KBspecR_matrix.RData"
  load(file.path(basepath,"RData",innameK))
  kbgood = which(KBspecR[,"opHi"]>0);  nkbgd = length(kbgood)
  # for (j in 1:nkbgd){print(specR[[KBspecR[kbgood[j],1]]]$strainOps[KBspecR[kbgood[j],2]:KBspecR[kbgood[j],3]])}
  
  ord.abKB = order(AbundKB$Average,decreasing=TRUE)
  # Relation between KB and AbundKB is derived with the following:-
  # Define jord.abund as an index on kbgood such that KB[kbgood[jord.abund]],] gives
  # the KB rows that have strains ordered in descending relative abundances. 
  # Derivation of jord.abund has been done maually.
  
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
  
  
  ## For the kth most abundant strain in KB that has strain data in specR here are the
  ## strain operon names:-
  #k = 5;    specRinds3 =  KBspecR[kbgood[jord.abund[k,1]],]
  #specR[[specRinds3[1]]]$strainOps[specRinds3[2]:specRinds3[3]]
  ## and here is the full specR detail:-
  #print(specR[[specRinds3[1]]])
  ra1 = AbundKB[ord.abKB[jord.abund[,2]],16]$Average
  RAkbg = ra1/sum(ra1)
  
  #print(cbind(KB[kbgood[jord.abund[,1]],8], AbundKB[ord.abKB[jord.abund[,2]],c(3,16)]))
  out = list(kbgood,ord.abKB,jord.abund,KB,specR,KBspecR,RAkbg)
}
######################################################################################
############################     RUN INITIALISATION      ##############################
#######################################################################################
subSampleFactor = 1
fastqPath = file.path(basepath,"fastq",whichSubunit)

constructing = FALSE
setwd(basepath)


#######################################################################################
#######################################################################################
##########################                               ##############################
##########################           MAIN BODY           ##############################
##########################                               ##############################
#######################################################################################
#######################################################################################
# PART 1: Specify the mock microbiome, and generate the fasta files of each operon in that
#         mock microbiome.
# Start with GFKB strain names of strains being used. These are manually revised based
# on specR operon names that have been trimmed to accord with blastn database requirements.

nmock = ifelse(whichMock=="mock50",50,59)
# DKB below gives  list(kbgood,ord.abKB,jord.abund,KB,specR,RAkbg)
DKB = linKBtospecR("specR_etal_grondDB.RData",constructing)
kbgood = DKB[[1]];              nkbgd = length(kbgood)
ord.abKB = DKB[[2]] 
jord.abund = DKB[[3]];          njab = nrow(jord.abund)
KB = DKB[[4]]
specR = DKB[[5]]
KBspecR = DKB[[6]]
RAkbg = DKB[[7]]
RA = RAkbg

# Derive the strain and species names of the organisms in the mock microbiome.
designStrainNames = sapply(1:njab,function(j){
  specRinds3 = KBspecR[kbgood[jord.abund[j,1]],]
  t0 = specR[[specRinds3[1]]]$strainOps[specRinds3[2]];  
  t1 = unlist(strsplit(t0,split=" "))
  out = paste(t1[1:(length(t1)-1)],sep=" ",collapse=" ")})

designSpeciesNames = sapply(1:nmock,function(j){
  specRinds3 = KBspecR[kbgood[jord.abund[j,1]],]
  out = specR[[specRinds3[1]]]$species.ordered })

# Generate matrix indOp giving the operon indices of the first and last operon of each strain.  
ind.strainOps = KBspecR[kbgood[jord.abund[,1]],"opLow"]
indupp.strainOps = KBspecR[kbgood[jord.abund[,1]],"opHi"]
num_ops_per_strain = sapply(1:length(ind.strainOps),function(j){indupp.strainOps[j] - ind.strainOps[j]+1})
maxopnum = sum(num_ops_per_strain)
indOp = matrix(c(c(1,1+cumsum(num_ops_per_strain[1:(nmock-1)])),cumsum(num_ops_per_strain)),nrow=nmock,ncol=2)

# Determine indices into specR of the speciesNames.
iSpec = KBspecR[kbgood[jord.abund[,1]]]
# Manually inspect the strainNames above and amend to align with those in specR.
# When complete, determine the index set for each strain within strainOps of the relevant specR element.
opIndex = vector("list",nmock)
for (j in 1:nmock){opIndex[[j]] = ind.strainOps[j]:indupp.strainOps[j] }
nOpsMock = sum(sapply(1:nmock,function(j){length(opIndex[[j]])}))
# Need to index from specR ops to GROND strain DB indexing.
jcount = 1
opSpecDBindex = matrix(0,nrow=length(specR),ncol=2)
for (j in 1:length(specR)){
  opSpecDBindex[j,1] = jcount
  nops = length(specR[[j]]$strainOps)
  jcount = jcount+nops
  opSpecDBindex[j,2] = jcount-1
}
# Now select the fasta sequences from the relevant multifasta<rRNAgene> sub-directory and write
# to outPath file.path(basepath,"fasta").
# Also write as text the sequence lengths of all the fasta files written. This is later used
# by the nanopore sequencer simulator, badread.
inPath = file.path(basepath,"GROND",paste("multifasta",whichSubunit,sep=""))
instem = paste("grondRefseq_",whichSubunit,"_op",sep="")
outPath = file.path(basepath,"fasta")
outstem = paste(whichMock,"_",whichSubunit,"_Op_",sep="")
fastaseqLengths = rep(0,10*nmock)


opcount = 0
OpIndex1=rep(0,nmock)
DBindex1=rep(0,nmock)
StrainOpSeqs = vector("list",nmock)
for (jsp in 1:nmock){
  jopSet = opSpecDBindex[iSpec[jsp],1]-1+opIndex[[jsp]]
  for (jjop in 1:length(jopSet)){
    jop = jopSet[jjop]
    opcount = opcount + 1
    cat("Index ", jsp,"Species ",designStrainNames[jsp]," strainDB op index ",jop,"  outname ",paste(outstem,opcount,".fasta",sep=""),"\n")
    if (jjop==1){
      OpIndex1[jsp] = opcount;    DBindex1[jsp] = jop
    }
    f1 = readFasta(dirPath=inPath,pattern=paste(whichSubunit,"_op",jop,sep="",collapse=""))
    # Modify ShortRead object's id to have no spaces, using "_" instead. Needed to ensure badread-generated
    # fastq objects have full mock strain operon name in header.
    newid = paste(unlist(strsplit(as.character(ShortRead::id(f1[1])),split=" ")),sep="_",collapse="_")
    newf1 = ShortRead(sread(f1[1]),id = BStringSet(newid))
    f1 = newf1
    if (width(f1) >5600) cat("\n Case jsp = ",jsp,"    jop = ",jop,"  too long.\n")
    writeFasta(f1,file=file.path(outPath,paste(outstem,opcount,".fasta",sep="")))
    fastaseqLengths[opcount] = width(f1)
  }
  cat("Num_ops_per_strain",num_ops_per_strain[jsp],"\n\n")
}
cat("Total operons: ",sum(num_ops_per_strain),"  Number of strains: ",nmock,"\n")
MockKB_OpsTable = data.frame(designStrainNames=designStrainNames,NopsperStrain=num_ops_per_strain, 
                                  OpIndex1=OpIndex1, DBindex1 = DBindex1)
print(MockKB_OpsTable)
outname = "MockKB_OpsTable_identify_strainOps_fastaSequences.RData"
save(MockKB_OpsTable,opSpecDBindex,file=file.path(basepath,"RData",outname))

fastaseqLengths = fastaseqLengths[1:opcount]
seqLenName = paste("frag",whichSubunit,"length.sh",sep="")
line1 = paste("lengths=(",paste(fastaseqLengths,sep=" ",collapse=" "),")",sep="")
writeLines(line1,con=file.path(basepath,"text",seqLenName))


# PART 2:  Use a set of system2() calls to execute badread with mainly fixed parameters, but allowing 
#          quality, operon number, and number of reads to be modified.  Also, note that
#          the name of the file of operon lengths may need to be varied in the shellScript being called.
# If the RA range is large - e.g. > 100 - it is very wasteful to generate the same number of reads per operon
# for all strains.  The following code reduces such wastage by adapting the number of reads generated per operon 
# using the planned relative abundances _assumed ordered largest to smallest - while ensuring that more
# than enough are generated.

RA.range = max(RA)/min(RA)
allTotfrags = 100 + 100*(maxTotfrags*1.5*RA%%100);    allTotfrags[1] = maxTotfrags
if (RA.range >10){
  cat("Large variation in relative abundances so use runs with different numbers of generated reads according to RA.\n")
  iord = order(RA,decreasing=TRUE)
  Ngen = round(maxTotfrags*1.2*RA[iord]/(RA[iord[1]]*100))*100 + 100;   Ngen[1] = maxTotfrags
  iset1 = which(Ngen>5000);   iset2 = setdiff(which(Ngen>1000),iset1); iset3 = setdiff(1:nmock,union(iset1,iset2))
  # Run each strain indexed by iset1 in its own slurm job.
  shellScriptName = "badread_mKB.sh"
  for (j1 in iset1){
    OpLow = indOp[j1,1];    OpHi = indOp[j1,2]
    total_frags = Ngen[j1]
    maxmemStr = "--mem=32GB";    maxtime = 1 + (OpHi - OpLow)*total_frags/150000 + 2
    maxtimeStr = paste("--time=",round(maxtime),":00:00",sep="")
    argstr1 = paste(maxmemStr, maxtimeStr, shellScriptName, basepath, whichSubunit,
                    OpLow, OpHi, identMean, identSD, total_frags,sep = " ")
    cat(argstr1,"\n")
    system2("sbatch",args=argstr1)
  }    #    end     j1    loop
  # Now set up the iset2 generation.
  OpLow = indOp[min(iset2),1];    OpHi = indOp[max(iset2),2]
  total_frags = Ngen[min(iset2)]
  maxmemStr = "--mem=32GB";    maxtime = 1 + (OpHi - OpLow)*total_frags/150000 + 2
  maxtimeStr = paste("--time=",round(maxtime),":00:00",sep="")
  argstr2 = paste(maxmemStr, maxtimeStr, shellScriptName, basepath, whichSubunit,
                  OpLow, OpHi, identMean, identSD, total_frags,sep = " ")
  cat(argstr2,"\n")
  system2("sbatch",args=argstr2)
  
  # And the iset3 generation - same pattern as iset2
  OpLow = indOp[min(iset3),1];    OpHi = indOp[max(iset3),2]
  total_frags = Ngen[min(iset3)]
  maxmemStr = "--mem=32GB";    maxtime = 1 + (OpHi - OpLow)*total_frags/150000 + 2
  maxtimeStr = paste("--time=",round(maxtime),":00:00",sep="")
  argstr3 = paste(maxmemStr, maxtimeStr, shellScriptName, basepath, whichSubunit,
                  OpLow, OpHi, identMean, identSD, total_frags,sep = " ")
  cat(argstr3,"\n")
  system2("sbatch",args=argstr3)
} else {   
  shellScriptName = "badread_mKB.sh"

  # Generate equal number of reads for each operon.
  # Generate Ns operon ranges of about equal size, then submit a slurm batch job for each.
  Ns = 6
  binSize = (nmock+round(Ns/2))%/%Ns
  for (k in 1:Ns){
    if (k<Ns){
      OpLow = binSize*(k-1) + 1;   OpHi = OpLow + binSize - 1
    } else {
      OpLow = OpHi+1;    OpHi = nmock
    }
    total_frags = maxTotfrags
    maxmemStr = "--mem=32GB";    maxtime = 1 + (indOp[OpHi,2] - indOp[OpLow,1])*total_frags/150000 + 2
    maxtimeStr = paste("--time=",round(maxtime),":00:00",sep="")
    argstr3 = paste(maxmemStr, maxtimeStr, shellScriptName, basepath, whichSubunit,
                    indOp[OpLow,1], indOp[OpHi,2], identMean, identSD, total_frags,sep = " ")
    cat(" Call:   sbatch",argstr3,"\n")
    system2("sbatch",args=argstr3)
  }
}

cat("PART 2 completed. \n")




