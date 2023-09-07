#!/usr/bin/env python

version = "v0.1.0.20220721"

import os,os.path,sys,re
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '/home/hongxiaoning/pipeline/bin'))
import myc
import time

def Parser_HighQuality(input_bed, output_bed):
    output = open(output_bed, 'w')
    with open(input_bed, 'r') as raw_bed:
        for line in raw_bed:
            ls = re.split("\t", line.strip())
            Discordants= int(ls[3])
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
            if  SplitReads  >= 2  \
            and Discordants >= 0 \
            and Score    >= 200 \
            and MC       >SD \
            and CC       <= 0.9:
                output.write(line)
            else:
                continue
    output.close()
            

if __name__ == '__main__':
    #Parameters to be input.
    parser=ArgumentParser(description="Detecting Circle DNA module, version {}".format(version))
    parser.add_argument("-s","--sample_bam", action="store", dest="sample_bam", 
            help="sample bam", required=True)
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
    tmp_dir = args.output_path +"/tmp"
    myc.createDirectory(tmp_dir)

    #get config information to dict
    config = myc.resolveConfig(args.cfg)
    prefix = args.prefix
    output_path = os.path.abspath(args.output_path)
    dirprefix = "{}/{}".format(output_path, prefix)

    productVersion = True
    if args.debug == True:
        productVersion = False
    # software and args
    samtools = config['samtools']
    cpu_circle = config["cpu_circle"]
    mem_circle = int(config["mem_circle"])/1000
    mem_circle = str(int(mem_circle)) + "G"

    # result note
    result_note_list = []

    # Preparing the files for Circle-Map
    #sample_bam = args.sample_bam
    sample_bam  = output_path + '/chr.{}.sorted.bam'.format(prefix)
    qname_bam  = output_path + '/qname.{}.sorted.bam'.format(prefix)
    #circular_bam=output_path + '/circular.{}.bam'.format(prefix)
    circular_sort_bam=output_path + '/circular.{}.sorted.bam'.format(prefix)
    circle_bed = output_path + '/{}.circle.bed'.format(prefix)
    circle_HighQuality_bed = output_path + '/{}.circle.HighQuality.bed'.format(prefix)
    CMD_script = output_path + '/CM_{}.sh'.format(prefix)

    result_note_list.append("Start time %s"%time.strftime("%Y-%m-%d %X", time.localtime()))
    result_note_list.append("START Circle-Map ")
    if productVersion or not os.path.exists(circle_bed):
        exportTMP_cmd= "export TMPDIR=" + tmp_dir #+ " && export PWD=" + output_path
        result_note_list.append(exportTMP_cmd)
        myc.run_command(exportTMP_cmd, ">>>Export TMPDIR Failed")
        Realign_cmd = " ".join([config["Circle-Map"], "Realign -t", cpu_circle, "-i", circular_sort_bam, "-qbam", qname_bam, "-sbam", sample_bam, \
                "-fasta", config['Ref'], "-o", circle_bed])#, "-dir", output_path])
        result_note_list.append(Realign_cmd)
        QC_cmd = " ".join([config["python"] ,config["Circle-Map_qc"], "-i", circle_bed, "-o", circle_HighQuality_bed])

        myc.write_readme(output_path + "/Circle-Map.README", version, result_note_list)
        myc.run_command(Realign_cmd, ">>>START Circle-Map Realign to detect the circular DNA Failed")
        Parser_HighQuality(circle_bed, circle_HighQuality_bed)
        result_note_list.append("Filter LowQuality circular DNA result")
        result_note_list.append(QC_cmd)
        #myc.run_command("QS -t "+cpu_circle+" -m "+mem_circle+" "+CMD_script, ">>>Qsub Failed")
        	
    # readme
    result_note_list.append("End time %s"%time.strftime("%Y-%m-%d %X", time.localtime()))
    myc.write_readme(output_path + "/Circle-Map.README", version, result_note_list)

    # remove files
    if productVersion:
        myc.rm_files("{}/*temp*".format(output_path))
        myc.rm_files("{}/*tmp*".format(output_path))
        myc.rm_files("{}/*bam*".format(output_path))

