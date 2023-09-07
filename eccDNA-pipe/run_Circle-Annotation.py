#!/usr/bin/env python

version = "v0.1.0.20220721"

from asyncore import write
from numbers import Rational
import os,os.path,sys,re,gzip,subprocess
from collections import defaultdict
from argparse import ArgumentParser
from tkinter.messagebox import NO
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '/home/hongxiaoning/pipeline/bin'))
import myc
import time,threading,multiprocessing
from concurrent.futures import ProcessPoolExecutor

def addtwodimdict(thedict, key_a, key_b, val):
    adict = thedict.keys()
    if key_a in adict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})

def Read_bed(input_bed):
    dict_stat = {}
    print ("Loading file: {} ...".format(input_bed))
    with open(input_bed, 'r') as raw_bed:
        for line in raw_bed:
            ls = re.split("\t", line.strip())
            Chrom = str(ls[0])
            Start = int(ls[1])
            End   = int(ls[2])
            Length= End - Start
            key = "_".join([Chrom, str(Start), str(End)])
            dict_stat[key] = [Chrom, Start, End, Length]
            SplitReads = int(ls[4])
            Score      = float(ls[5])
            MC         = float(ls[6]) #MC , Mean coverage
            SD         = float(ls[7]) #SD , Standard deviation
            try:
                CS     = float(ls[8]) #CS , Coverage increase in the start coordinate
            except ValueError:
                CS     = 0
            try:
                CE     = float(ls[9]) #CE , Coverage increase in the end coordinate
            except ValueError:
                CE     = 0
            try:
                CC     = float(ls[10])#CC , Coverage continuity
            except ValueError:
                CC     = 1
    return dict_stat
            
def Read_gff(gff_file):
    dict_gff = {}
    print ("Loading file: {} ...".format(gff_file))
    if gff_file.endswith(".gz") == False:
        f_gff = open(gff_file,'r')
    else:
        f_gff = gzip.open(gff_file,'rt')
    for line in f_gff:
        line=line.rstrip('\n')
        if line and (not line.startswith('#')):
            lines=line.split('\t')
            if lines[2]=='gene':
                gene =re.search('gene_name=([^;]+)' ,lines[8]).group(1) 
                chrom=lines[0]
                start=lines[3]
                end  =lines[4] 
                key  ="_".join([chrom, start, end, gene]) 
                dict_gff[key] = [chrom, int(start), int(end), gene]
    
    return dict_gff
            
def Read_miRNA(gff_file):
    dict_gff = {}
    print ("Loading file: {} ...".format(gff_file))
    if gff_file.endswith(".gz") == False:
        f_gff = open(gff_file,'r')
    else:
        f_gff = gzip.open(gff_file,'rt')
    for line in f_gff:
        line=line.rstrip('\n')
        if line and (not line.startswith('#')):
            lines=line.split('\t')
            if lines[2]=='miRNA':
                gene =re.search('Name=([^;]+)' ,lines[8]).group(1) 
                chrom=lines[0]
                start=lines[3]
                end  =lines[4] 
                key  ="_".join([chrom, start, end, gene]) 
                dict_gff[key] = [chrom, int(start), int(end), gene]
    
    return dict_gff

def Read_cCREs(input_bed):
    dict_enhancer = {}
    print ("Loading file: {} ...".format(input_bed))
    with open(input_bed, 'r') as raw_bed:
        for line in raw_bed:
            ls = re.split("\t", line.strip())
            Chrom = str(ls[0])
            Start = int(ls[1])
            End   = int(ls[2])
            Length= End - Start
            key = str(ls[3] )
            dict_enhancer[key] = [Chrom, Start, End, key]
    return dict_enhancer

def Read_TE(input_bed):
    dict_TE = {}
    print ("Loading file: {} ...".format(input_bed))
    with open(input_bed, 'r') as raw_bed:
        for line in raw_bed:
            ls = re.split("\t", line.strip())
            Chrom = str(ls[0])
            Start = int(ls[1])
            End   = int(ls[2])
            key = ls[4]+":"+ls[0]+":"+ls[1]+"-"+ls[2]
            dict_TE[key] = [Chrom, Start, End, key]
    return dict_TE

def ann_gene(gene_gff,bed_file,output_path):
    gff = Read_gff(gene_gff)
    sample_info = Read_bed(bed_file)
    sample= re.split("\.",re.split("\/",bed_file)[-1])[0]
    dict_rmdup_eccDNA = {}
    sample_gene_output = open(output_path+"/"+sample + ".Gene.txt", 'w')
    sample_gene_output.write("Sample\tChrom\tStart\tEnd\tGene_Name\tGene_Start\tGene_End\tOverlap_Length\n")
    for key in sorted(sample_info):
        Chrom = sample_info[key][0]
        Start = sample_info[key][1]
        End   = sample_info[key][2]
        Length= sample_info[key][3]
        #Gene Annotation
        for value in gff.values():
            gene_chr  = value[0]
            gene_start= value[1]
            gene_end  = value[2]
            gene_name = value[3]
            if Chrom == gene_chr:
                overlap_length = max(0, min(End, gene_end) - max(gene_start, Start) +1)
                if overlap_length >=70:
                    sample_gene_output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(sample,Chrom,Start,End,gene_name,gene_start,gene_end,overlap_length))
                    dict_rmdup_eccDNA[key] = key
            else:
                continue
    sample_gene_output.close()
    with open(output_path+"/"+"result.stat", 'a') as f_output2:
        eccDNA_Gene_number = len(dict_rmdup_eccDNA)
        f_output2.write("%s\tGene\t%s\n" %(sample, eccDNA_Gene_number))        

def ann_miRNA(miRNA_gff,bed_file,output_path):
    gff = Read_miRNA(miRNA_gff)
    sample_info = Read_bed(bed_file)
    sample= re.split("\.",re.split("\/",bed_file)[-1])[0]
    dict_rmdup_miRNA = {}
    sample_miRNA_output = open(output_path+"/"+sample + ".miRNA.txt", 'w')
    sample_miRNA_output.write("Sample\tChrom\tStart\tEnd\tmiRNA_Name\tmiRNA_Start\tmiRNA_End\tOverlap_Length\n")
    for key in sorted(sample_info):
        Chrom = sample_info[key][0]
        Start = sample_info[key][1]
        End   = sample_info[key][2]
        Length= sample_info[key][3]
        #miRNA Annotation
        for value in gff.values():
            gene_chr  = value[0]
            gene_start= value[1]
            gene_end  = value[2]
            gene_name = value[3]
            if Chrom == gene_chr:
                overlap_length = max(0, min(End, gene_end) - max(gene_start, Start) +1)
                if overlap_length >= 15:
                    sample_miRNA_output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(sample,Chrom,Start,End,gene_name,gene_start,gene_end,overlap_length))
                    dict_rmdup_miRNA[key] = key
            else:
                continue
    sample_miRNA_output.close()
    with open(output_path+"/"+"result.stat", 'a') as f_output2:
        eccDNA_miRNA_number = len(dict_rmdup_miRNA)
        f_output2.write("%s\tmiRNA\t%s\n" %(sample, eccDNA_miRNA_number))        

def ann_cCREs(cCREs_bed,bed_file,output_path):
    #gff = Read_cCREs(cCREs_bed)
    #sample_info = Read_bed(bed_file)
    sample= re.split("\.",re.split("\/",bed_file)[-1])[0]
    dict_rmdup_cCREs = {}
    sample_gene_output = open(output_path+"/"+sample + ".cCREs.txt", 'w')
    sample_gene_output.write("Sample\tChrom\tStart\tEnd\tcCRE_ID\tcCRE_Start\tcCRE_End\tOverlap_Length\n")
    cmd = " ".join(["bedtools intersect -a", bed_file, "-b", cCREs_bed, "-wo"])
    returnvalue,output = subprocess.getstatusoutput(cmd)
    list_res=re.split("\n", output)
    for line  in list_res:
        ls = re.split("\t",line)
        overlap_length = int(ls[-1])+1
        if overlap_length >= 70:
            Chrom = ls[0]
            Start = ls[1]
            End   = ls[2]
            eccDNA_id = Chrom+":"+str(Start)+"-"+str(End)
            cCRE_start= ls[-6]
            cCRE_end  = ls[-5]
            cCRE_name=Chrom+":"+str(cCRE_start)+"-"+str(cCRE_end)
            sample_gene_output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(sample,Chrom,Start,End,cCRE_name,cCRE_start,cCRE_end,overlap_length))
            dict_rmdup_cCREs[eccDNA_id] = 1

    sample_gene_output.close()
    with open(output_path+"/"+"result.stat", 'a') as f_output2:
        cCRE_number = len(dict_rmdup_cCREs)
        f_output2.write("%s\tcCRE\t%s\n" %(sample, cCRE_number))        

def ann_TE(TEs_bed,bed_file,output_path):
    #TE = Read_TE(TEs_bed)
    #sample_info = Read_bed(bed_file)
    sample= re.split("\.",re.split("\/",bed_file)[-1])[0]
    
    dict_rmdup_TE = {}
    sample_TE_output = open(output_path+"/"+sample + ".TEs.txt", 'w')
    sample_TE_output.write("Sample\tChrom\tStart\tEnd\tTE_ID\tTE_Start\tTE_End\tOverlap_Length\n")
    # Annotation
    cmd = " ".join(["bedtools intersect -a", bed_file, "-b", TEs_bed, "-wo"])
    returnvalue,output = subprocess.getstatusoutput(cmd)
    list_res=re.split("\n", output)
    for line  in list_res:
        ls = re.split("\t",line)
        overlap_length = int(ls[-1])+1
        if overlap_length >= 70:
            Chrom = ls[0]
            Start = ls[1]
            End   = ls[2]
            eccDNA_id = Chrom+":"+str(Start)+"-"+str(End)
            TE_start= ls[-5]
            TE_end  = ls[-4]
            TE_name=Chrom+":"+str(TE_start)+"-"+str(TE_end)
            sample_TE_output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(sample,Chrom,Start,End,TE_name,TE_start,TE_end,overlap_length))
            dict_rmdup_TE[eccDNA_id] = 1

        
    sample_TE_output.close()
    with open(output_path+"/"+"result.stat", 'a') as f_output2:
        eccDNA_TE_number = len(dict_rmdup_TE)
        f_output2.write("%s\tTE\t%s\n" %(sample, eccDNA_TE_number))

def dataset_coverge(dataset_bed, bed_file):
    cmd = "bedtools coverage  -a " + dataset_bed + " -b " + bed_file + "|awk '{sum+=$NF} END {print sum/NR}'"
    returnvalue,output = subprocess.getstatusoutput(cmd)
    
    return(float(output))

def dataset_cover_number(dataset_bed, bed_file):
    cmd = "bedtools coverage  -b " + dataset_bed + " -a " + bed_file 
    returnvalue,output = subprocess.getstatusoutput(cmd)
    number = 0
    list_res=re.split("\n", output)
    for line  in list_res:
        line = re.split("\t",line)
        if int(line[-3]) > 70:
            number +=1

    return(number)

def Normalized_dataset_coverage(Chr_bed, Type_bed, sample_bed):
    sample_ratio = dataset_coverge(Type_bed,sample_bed)
    Type_ratio = dataset_coverge(Chr_bed,Type_bed)
    Normalized_ratio=round(sample_ratio/Type_ratio, 4)
    
    return(Normalized_ratio)

def Chromosome_eccDNA_Number(sample_bed, config, output_path):
    cmd = "bedtools coverage  -a " + config['Chr_bed'] + " -b " + sample_bed 
    returnvalue,output = subprocess.getstatusoutput(cmd)
    list_res=re.split("\n", output) 
    sample= re.split("\.",re.split("\/",sample_bed)[-1])[0]
    f_output = open(output_path+"/Chromosome_Number.stat", 'a')
    for line  in list_res:
        line = re.split("\t",line)
        eccDNA_number = int(line[-4])
        Chromosome_size= int(line[-2])
        number = round(eccDNA_number/(Chromosome_size/1000000),4) #per Mb
    
        f_output.write("%s\t%s\t%s\n" %(sample, line[0], number))        
    f_output.close()

def ann_genomic_coverage(sample_bed,config,output_path):
    sample= re.split("\.",re.split("\/",sample_bed)[-1])[0]
    with open(output_path+"/Genomic_coverage.stat", 'a') as f_output:
        Normalized_3UTR_ratio=Normalized_dataset_coverage(config["chr_bed"], config["3UTR_bed"], sample_bed)
        f_output.write("%s\t3UTR\t%s\n" %(sample, Normalized_3UTR_ratio))        
        Normalized_5UTR_ratio=Normalized_dataset_coverage(config["chr_bed"], config["5UTR_bed"], sample_bed)
        f_output.write("%s\t5UTR\t%s\n" %(sample, Normalized_5UTR_ratio))        
        Normalized_CpG_ratio=Normalized_dataset_coverage(config["chr_bed"], config["CpGIsland_bed"], sample_bed)
        f_output.write("%s\tCpG\t%s\n" %(sample, Normalized_CpG_ratio))        
        Normalized_Exon_ratio=Normalized_dataset_coverage(config["chr_bed"], config["Exon_bed"], sample_bed)
        f_output.write("%s\tExon\t%s\n" %(sample, Normalized_Exon_ratio))        
        Normalized_Gene2KbD_ratio=Normalized_dataset_coverage(config["chr_bed"], config["Gene2KbD_bed"], sample_bed)
        f_output.write("%s\tGene2KD\t%s\n" %(sample, Normalized_Gene2KbD_ratio))        
        Normalized_Gene2KbU_ratio=Normalized_dataset_coverage(config["chr_bed"], config["Gene2KbU_bed"], sample_bed)
        f_output.write("%s\tGene2KU\t%s\n" %(sample, Normalized_Gene2KbU_ratio))        
        Normalized_Intron_ratio=Normalized_dataset_coverage(config["chr_bed"], config["Introns_bed"], sample_bed)
        f_output.write("%s\tIntron\t%s\n" %(sample, Normalized_Intron_ratio))        

if __name__ == '__main__':
    #Parameters to be input.
    parser=ArgumentParser(description="Circle DNA Annotation module, version {}".format(version))
    parser.add_argument("-s","--sample_bed", action="store", dest="sample_bed", 
            help="sample bed", required=True)
    parser.add_argument("-o","--output_path", action="store", dest="output_path", 
            help="output path", required=True)
    parser.add_argument("-c","--cfg", action="store", dest="cfg", 
            help="pipeline configure file", required=True)
    parser.add_argument("-p","--prefix", action="store", dest="prefix", 
            help="output prefix, sample id", required=True)
    parser.add_argument("--debug", action="store_true", dest="debug", 
            help="for debug module, skipping existed file, default False", default=False)

    args = parser.parse_args()

    #check and create dir
    myc.check_files([args.cfg])
    myc.createDirectory(args.output_path)
    #end check and create dir

    #get config information to dict
    config = myc.resolveConfig(args.cfg)
    prefix = args.prefix
    output_path = os.path.abspath(args.output_path)
    dirprefix = "{}/{}".format(output_path, prefix)

    productVersion = True
    if args.debug == True:
        productVersion = False

    # result note
    result_note_list = []

    # Preparing the files for Circle-Map
    sample_bed = args.sample_bed
    sample_eccDNA= output_path+"/{}.eccDNA.txt".format(prefix)
    
    if productVersion or not os.path.exists(sample_eccDNA):
        output_eccDNA = open(sample_eccDNA, 'w')
        output_eccDNA.write("Sample\tChrom\tStart\tEnd\tLength\n")
        dict_total_bed = defaultdict(int)
        
        result_note_list.append("Start time %s"%time.strftime("%Y-%m-%d %X", time.localtime()))
        result_note_list.append("Run Circle-Annotation ")
        #read bed file
        sample_info = Read_bed(sample_bed)
        sample= re.split("\.",re.split("\/",sample_bed)[-1])[0]
        dict_total_bed[sample] = len(sample_info.keys())
        for key in sorted(sample_info):
            Chrom = sample_info[key][0]
            Start = sample_info[key][1]
            End   = sample_info[key][2]
            Length= sample_info[key][3]
            output_eccDNA.write("%s\t%s\t%s\t%s\t%s\n" %(sample, Chrom, Start, End, Length))        

        output_eccDNA.close()
        with open(output_path+"/result.stat", 'a') as output_stat:
            for sample in sorted(dict_total_bed):
                output_stat.write("%s\teccDNA\t%s\n" %(sample, dict_total_bed[sample]))        
        
        #Chrmosome eccDNA Number Per Mb 
        Chromosome_eccDNA_Number(sample_bed, config, output_path)
        #Normalized Genomic Coverage 
        ann_genomic_coverage(sample_bed, config, output_path)
        #gene 
        T1=multiprocessing.Process(target=ann_gene,args=(config['gene_gff'],sample_bed,output_path,))
        #miRNA
        T2=multiprocessing.Process(target=ann_miRNA,args=(config['miRNA_gff'],sample_bed,output_path,))
        #TE
        T3=multiprocessing.Process(target=ann_TE,args=(config['TEs_bed'],sample_bed,output_path,))
        #cCREs
        T4=multiprocessing.Process(target=ann_cCREs,args=(config['cCREs_bed'],sample_bed,output_path,))

        T1.start()
        T2.start()
        T3.start()
        T4.start()


        	
    # readme
    result_note_list.append("End time %s"%time.strftime("%Y-%m-%d %X", time.localtime()))
    myc.write_readme(output_path + "/Circle-Annotation-%s.README"%prefix, version, result_note_list)

    # remove files
    if productVersion:
        myc.rm_files("{}/*temp*".format(output_path))


