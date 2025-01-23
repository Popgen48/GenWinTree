import re
import os
import logging
from pysam import VariantFile
from collections import OrderedDict
import glob
import json

# Load the configuration file
with open("config.json") as config_file:
    config = json.load(config_file)

# Paths and Parameters
RaxML= config["tools"]["RaxML"]
inDir= config["base_dir"]["inputDir"]
outDir= config["base_dir"]["outputDir"]
windowSize= int(config["params"]["windowSize"])
stepSize= int(config["params"]["stepSize"])
threads= int(config["params"]["threads"])
sampList= config["params"]["sampleList"]
vcfPattern= config["params"]["vcfPattern"]
genoDict= {
    "AG": "R", "GA": "R", "CT": "Y", "TC": "Y", "GC": "S", "CG": "S",
    "AT": "W", "TA": "W", "GT": "K", "TG": "K", "AC": "M", "CA": "M",
    "AA": "A", "GG": "G", "CC": "C", "TT": "T"
}


# Ensure the output directory exists
os.makedirs(outDir, exist_ok=True)

# Logging setup
LOG_FILE = os.path.join(outDir, "workflow.log")
logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logging.info("Starting phylogenetic workflow.")

# Extract chromosome names
CHROM, = glob_wildcards(inDir+"{chrom}"+vcfPattern)


# Snakemake Rules
rule all:
    input:
        expand(outDir+"{chrom}"+"{chrom}.tree.out",chrom=CHROM)

rule makeTrees:
    input:
        vcfIn=inDir+"{chrom}"+vcfPattern
    output:
        outDir+"{chrom}"+"{chrom}.tree.out"
    params:
        fastaIn=outDir+"{chrom}"+"{chrom}.fasta",
        treeOut=outDir+"{chrom}"+"{chrom}.fasta.raxml.bestTree"
    threads:
        threads
    run:
        chrmName=wildcards.chrom[1:]
        vcfInFile=VariantFile(input.vcfIn)
        for i,v in enumerate(list(vcfInFile.header.records)):
                if "##contig" in str(v):
                    pattern=re.compile(r"##contig\=\<ID\=([^,]+)\,length\=([0-9]+)\,(.*)")
                    match=re.findall(pattern,str(v))
                    if match[0][0]==chrmName:
                        chrmLength=int(match[0][1])
        for i in range(1,chrmLength,stepSize):
            j=i+windowSize
            if i+windowSize>=chrmLength:
                j=chrmLength
            tmpSeqDict=OrderedDict()
            totalCount=1
            tmpGenotypedDict=OrderedDict()
            tmpDepthDict=OrderedDict()
            for sampleName in sampList:
                tmpSeqDict[sampleName]=""
                tmpGenotypedDict[sampleName]=1
                tmpDepthDict[sampleName]=0
            for rec in vcfInFile.fetch(chrmName,i,j):
                totalCount+=1
                tmpGenoDict={0:rec.ref,1:rec.alts[0]}
                for sampleName in sampList:
                    if rec.samples[sampleName]["GT"][0]==None:
                        tmpSeqDict[sampleName]+="N"
                    else:
                        tmpGenotypedDict[sampleName]+=1
                        tmpDepthDict[sampleName]+=rec.samples[sampleName]["DP"]
                        tmpGeno=tmpGenoDict[rec.samples[sampleName]["GT"][0]]+tmpGenoDict[rec.samples[sampleName]["GT"][1]]
                        tmpSeqDict[sampleName]+=genoDict[tmpGeno]
            outP=str(params.fastaIn)
            destF=open(outP,"w")
            for sampleName in sampList:
                destF.write(">"+sampleName)
                destF.write("\n")
                destF.write(tmpSeqDict[sampleName])
                destF.write("\n")
            destF.close()
            genotypedList=list(tmpGenotypedDict.values())
            depthList=list(tmpDepthDict.values())
            avgDepthList=[str(round(depthList[i]/genotypedList[i],3)) for i in range(len(genotypedList))]
            avgGenotypList=[round(genotypedList[i]/totalCount,3) for i in range(len(genotypedList))]
            if min(avgGenotypList)>=0.95 and totalCount>11:
                shell(RaxML+" --parse --msa {params.fastaIn} --model GTR+G --threads {threads}")
                shell(RaxML+" --msa {params.fastaIn} --model GTR+G --threads {threads} --out ERR454947urial --seed 2 --force perf_threads --redo")
                outRes=str(output)
                destRes=open(outRes,"a")
                outT=str(params.treeOut)
                localTree=""
                with open(outT) as source1:
                    for line in source1:
                        localTree=line.rstrip()
                avgGenotypList=map(str,avgGenotypList)
                destRes.write(chrmName+"\t"+str(i)+"\t"+str(j)+"\t"+str(totalCount)+"\t"+\
                ",".join(avgGenotypList)+"\t"+",".join(avgDepthList)+"\t"+localTree+"\n")
                destRes.close()
                shell("rm {params.fastaIn}")
		
        logging.info(f"Finished processing chromosome: {chrmName}")

