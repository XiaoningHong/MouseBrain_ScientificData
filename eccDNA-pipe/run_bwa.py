#!/usr/bin/env python

version = "v0.1.0.20220610"

import os,os.path,sys,re
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '/home/hongxiaoning/pipeline/bin'))
import myc
import time
import pysam

def ReadCount(InputFile, OutputFile):
    if not os.path.isfile(InputFile):
        print ('File does not exist at location provided, please check (tried "{}")'.format(InputFile))
        exit(2)
    print ("Loading file: {} ...".format(InputFile))  
    f_output = open(OutputFile, 'w')
    samfile = pysam.AlignmentFile(InputFile, "rb")
    result  = samfile.get_index_statistics()    
    for line in result:
        match_obj = re.search(r'contig=\'(.+)\', mapped=(\d+)', str(line))
        contig = match_obj.group(1)
        mapped = match_obj.group(2)
        f_output.write('%s\t%s\n' %(contig, mapped))
        
    f_output.close()
    

if __name__ == '__main__':
    #Parameters to be input.
    parser=ArgumentParser(description="RNA mapping module, version {}".format(version))
    parser.add_argument("-fqs","--fqs", action="store", dest="fqs", 
            help="fastq, for paired end, split by comma", required=True)
    parser.add_argument("-o","--output_path", action="store", dest="output_path", 
            help="output path", required=True)
    parser.add_argument("-c","--cfg", action="store", dest="cfg", 
            help="pipeline configure file", required=True)
    parser.add_argument("-p","--prefix", action="store", dest="prefix", 
            help="output prefix, sample id", required=True)
    parser.add_argument("--markDuplicate", action="store", dest="markDuplicate", 
            help="do markDuplicate, default no", default="no", choices = ["no", "markDup", "removeDup"])
    parser.add_argument("--debug", action="store_true", dest="debug", 
            help="for debug module, skipping existed file, default False", default=False)

    args = parser.parse_args()

    #check and create dir
    myc.check_files([args.cfg])
    myc.createDirectory(args.output_path)

    #get config information to dict
    config = myc.resolveConfig(args.cfg)
    prefix = args.prefix
    output_path = os.path.abspath(args.output_path)
    dirprefix = "{}/{}".format(output_path, prefix)

    productVersion = True
    if args.debug == True:
        productVersion = False

    # software and args
    java = config["java_latest"]
    picard = config["picard_latest"]
    samtools = config['samtools']
    cpu_aln = config["cpu_aln"]
    mem_aln = int(config["mem_aln"])/1000
    mem_aln = str(int(mem_aln)) + "G"
    #tmp_dir = config["tmp_dir"]

    # result note
    result_note_list = []

    fqs = args.fqs
    fq_list = fqs.split(",")
    myc.check_files(fq_list)

    # BWA-MEM Mapping
    sort_bam = "{}.sort.bam".format(dirprefix)
    sample_output = '{}.ReadCount.txt'.format(dirprefix)
    group_header = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:WGS".format(sample = prefix)

    result_note_list.append("Start time %s"%time.strftime("%Y-%m-%d %X", time.localtime()))
    result_note_list.append("STAR mapping and sort")
    
    BWA_cmd = " ".join([config['bwa'], "mem -M -t", cpu_aln, "-R", "\""+group_header+"\"", config['Ref'], fq_list[0], fq_list[1], "|samtools sort -@", cpu_aln, "-o", sort_bam, "-"])
    ReadCount(sort_bam, sample_output)
    
    if productVersion or not os.path.exists(bam):
        myc.rm_files("{}/*tmp*".format(output_path))
        result_note_list.append(BWA_cmd)
        myc.run_command(BWA_cmd, ">>>BWA maping Failed")
        sort_cmd = "set -e;{} index {};{} flagstat {} >{}.stat;set +e".format(samtools,sort_bam,samtools,sort_bam,sort_bam)
        result_note_list.append(sort_cmd)
        myc.run_command(sort_cmd, ">>>Samtools sort and index Failed")
    else:
        if not os.path.isfile(sample_output):
            myc.run_command(ReadCount_cmd, ">>>ReadCount Failed")

    # mark duplicate if needs
    if args.markDuplicate != "no":
        input_bam = sort_bam
        if args.markDuplicate == "removeDup":
            bam = "{}.dedup.bam".format(dirprefix)
            remove_flag = "REMOVE_DUPLICATES=true"
            result_note_list.append("Picard mark duplicate and remove duplicate")
        else:
            bam = "{}.dedup.bam".format(dirprefix)
            remove_flag = "REMOVE_DUPLICATES=false"
            result_note_list.append("Picard mark duplicate")

        if productVersion or not os.path.exists(bam):
            tmp_dir='{}.duptmp'.format(dirprefix)
            cmd_dup = " ".join([java, "-jar -Xmx" + mem_aln, picard, "MarkDuplicates", "TMP_DIR=" + tmp_dir,
                "I=" + input_bam, "O=" + bam, "CREATE_INDEX=true", "M={}.dedup.metrics".format(dirprefix),
                remove_flag])
            myc.run_command(cmd_dup, ">>>picard mark duplicate Failed")
            result_note_list.append(cmd_dup)

    # readme
    result_note_list.append("End time %s"%time.strftime("%Y-%m-%d %X", time.localtime()))
    myc.write_readme(output_path + "/aln.README", version, result_note_list)

    # remove files
    if productVersion:
        if args.markDuplicate != "no":            
            myc.rm_files("{}.sorted.bam*".format(dirprefix))
            myc.rm_files("{}.duptmp".format(dirprefix))
        myc.rm_files("{}*tmp*".format(dirprefix))

            
