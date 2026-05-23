using Pkg
Pkg.add("ArgParse")
Pkg.add(url="https://github.com/MurrellGroup/NextGenSeqUtils.jl");
Pkg.add(url="https://github.com/MurrellGroup/DPMeansClustering.jl");
Pkg.add(url="https://github.com/MurrellGroup/RobustAmpliconDenoising.jl");
Pkg.add(url="https://github.com/MurrellGroup/SeqUMAP.jl");
Pkg.add("DataFrames");
Pkg.add("DelimitedFiles");
Pkg.add("StatsBase")
Pkg.add("PyPlot");

#using ArgParse
#commandLine=ArgParseSettings()
#println(commandLine)
#@add_arg_table! commandLine begin
#    "basepath"
#        help = "The parent directory for all sub-directories required for the processing pipeline"
#        arg_type = AbstractString
#        required = true
#    "--whichMock", "-m"
#        help = "Name of the mock microbiome - e.g. mockKB, mSerF,mSerS2,mSriSA3,mSriSZ4"
#        arg_type = AbstractString
#        required = true
#    "--whichSubunit", "-s"
#        help = "16S, 23S, or rrn"
#        arg_type = AbstractString
#        required = true
#    "--whichCase", "-c"
#        help = "which set of quality parameters, e.g. 03, 10, 11"
#        arg_type = AbstractString
#        required = true
#    "whichSubMock"
#        help = "For Sereika data, indexes the sub-samples.  default = 0"
#        arg_type = Int
#        default = 0 
#    "whichPair"
#        help = "For Srinivas data specifying which primer pair"
#        arg_type = Int
#        default = 0
#   "uppThreshErrRate"
#        help = "Upper threshold on read error rate for pre-RAD quality filtering. Passed as an integer, being the reciprocal of the error rate."
#        arg_type = AbstractString
#        default = 0.01
#end

basepath = ARGS[1]
wMock = ARGS[2]
wSub = ARGS[3]
wCase = ARGS[4]
wSM = ARGS[5]
wP = ARGS[6]
upperERstring = ARGS[7]

#parsed_args = parse_args(commandLine)
#basepath=parsed_args["basepath"]
#wMock=parsed_args["whichMock"]
#wSub=parsed_args["whichSubunit"]
#wCase=parsed_args["whichCase"]
#wSM=parsed_args["whichSubMock"]
#wP=parsed_args["whichPair"]
#upperERstring = parsed_args["uppThreshErrRate"]
println(typeof(upperERstring))
uppER=parse(Float64,upperERstring)
println(uppER)
ERstr = SubString(string(uppER),1,minimum([6,length(string(uppER))]))
ERchar = collect("0.0000")
charER = collect(ERstr)
for j in 1:6
  ERchar[j] = (j>length(charER) ? "0"[1] : charER[j])
end
ERstr = join(ERchar,"")
ERstring = join(["ER_less_",ERstr],"")
    
println("Arguments table constructed.")
println("basepath:         ", basepath)
println("whichMock:        ", wMock)
println("whichSubunit:     ", wSub)
println("whichCase:        ", wCase )
println("whichSubMock:     ", wSM)
println("whichPair:        ", wP)
println("uppThreshErrRate: ",uppER)

println(length(wMock))
if wMock=="mockKB"
  mstr = "mockKB"
  dataset = join(["mKB",wSub,join(["C",wCase],"")],"_")
  stem1 = join([wMock,wSub,"Case",wCase],"_")
elseif wMock=="mSerF"
  mstr = "mSerF"
  dataset = join([wMock,wSub,join(["C",wCase],"")],"_")
  stem1 = join([wMock,wSub,"Case",wCase],"_")
elseif SubString(wMock,1,5) == "mSerS"
  mstr = "mSerS"
  dataset = join([wMock,wSub,join(["C",wCase],"")],"_")
  stem1 = join([wMock,wSub,"Case",wCase],"_")
elseif SubString(wMock,1,6)=="mSriSA"
  mstr = "SA"
  dataset = join([join([mstr,wP],""),wSub,join(["C",wCase],"")],"_")
  stem1 = join([join([mstr,wP],""),wSub,"Case",wCase],"_")
elseif SubString(wMock,1,6)=="mSriSZ"
  mstr = "SZ"
  dataset = join([join([mstr,wP],""),wSub,join(["C",wCase],"")],"_")
  stem1 = join([join([mstr,wP],""),wSub,"Case",wCase],"_")
else
  mstr = "invalid"
  dataset = "None"
  stem1 = "None"
end


println(dataset)
println(stem1)

stem2 = "filtered_denoise"
suffix1 = "out.txt"
fname1 = join([stem1,".fastq"],"")
fastq1 = joinpath(basepath,"fastq",fname1)
fname2 = join([join([stem1,ERstring,"filtered"],"_"),".fastq"],"")
fastq2 = joinpath(basepath,"fastq",fname2)

using NextGenSeqUtils,DPMeansClustering,RobustAmpliconDenoising

if wSub == "16S"
  LowLen = 1100;   HiLen = 1600
elseif wSub == "23S"
  LowLen = 2100;   HiLen = 2950
else
  LowLen = 4300;   HiLen = 5600 
end 

using PyPlot
# Length and quality filter the data and then denoise it. Write the templates (==ASVs) to a .fasta file.
fastq_filter(fastq1,fastq2,error_rate=uppER,min_length=LowLen,max_length=HiLen)
length_vs_qual(fastq2)
plotname1 = join([join(["errorRate_vs_length",dataset,ERstring],"_"),".pdf"],"")
println(plotname1)
savefig(joinpath(basepath,"plots",plotname1))
reads,phreds,names = read_fastq(fastq2);
templates,template_sizes,template_indices=denoise(reads);
if isempty(templates) 
  println("No ASVs formed - no further computation.")
else 
  begin
    fname3 = join([join([stem1,ERstring,"filtered"],"_"),".fasta"],"")
    println(fname3)
    write_fasta(joinpath(basepath,"fasta",fname3),templates,names=["seq_$(j)_$(template_sizes[j])" for j in 1:length(template_sizes)])

    # Save as text files the names of the retained reads, their indices in the input file,
    # and the UMAP coordinates of the ASVs and the retained reads.
    using DataFrames, DelimitedFiles
    fnameN = join([stem1,ERstring,stem2,"names",suffix1],"_")
    open(joinpath(basepath,"text",fnameN), "w") do file
        writedlm(file,names)
    end

    fnameI = join([stem1,ERstring, stem2, "indices", suffix1],"_")
    open(joinpath(basepath,"text",fnameI), "w") do file
        writedlm(file,template_indices)
    end

    fnameT = join([stem1, ERstring, stem2, "templates", suffix1],"_")
    open(joinpath(basepath,"text",fnameT), "w") do file
        writedlm(file,templates)
    end

    # Generate UMAP coordinates for the combination of ASVs and the reads from which they were derived.
    using SeqUMAP
    fnameP = join([stem1,ERstring, stem2, "proj", suffix1],"_")
    proj = sequmap(vcat(templates,reads),2, k=6, n_neighbors=10, pca=false, min_dist=1.0);
    open(joinpath(basepath,"text",fnameP), "w") do file
        writedlm(file,proj)
    end
    println("UMAP2D projection calculation completed.")

    fig, ax = subplots()
    ax.scatter(proj[1,(length(templates)+1):end],proj[2,(length(templates)+1):end]; color="y", s=2.0, linewidth = 0.0)
    ax.scatter(proj[1,1:length(templates)],proj[2,1:length(templates)]; marker="+", c="black", s=13.0, linewidth=1.0)
    plottitle = join(["UMAP",dataset,"+ve strand; reads(yellow), ASVs(blk,+): ErrRate<",round(uppER,digits=4)]," ");
    title(plottitle)
    fnameUM = join([join(["umap2D",dataset, ERstring],"_"),".pdf"],"");
    println(fnameUM)
    savefig(joinpath(basepath,"plots",fnameUM))
    println("FINISHED")
  end
end
