#!/usr/bin/env python

version = "v0.1.0.20220721"

import os,os.path,sys,re
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '/home/hongxiaoning/pipeline/bin'))
import myc
import time


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
    bedtools = config['bedtools']
    cpu_circle = str(4)
    #cpu_circle = config["cpu_circle"]
    mem_circle = int(config["mem_circle"])/1000
    mem_circle = str(int(mem_circle)) + "G"

    # result note
    result_note_list = []

    # Preparing the files for Circle-Map
    #sample_bam = args.sample_bam
    sample_bam = output_path + '/chr.{}.sorted.bam'.format(prefix)
    qname_bam  = output_path + '/qname.{}.sorted.bam'.format(prefix)
    circular_bam=output_path + '/circular.{}.bam'.format(prefix)
    circular_sort_bam=output_path + '/circular.{}.sorted.bam'.format(prefix)

    result_note_list.append("Start time %s"%time.strftime("%Y-%m-%d %X", time.localtime()))
    result_note_list.append("STAR Prepare for Circle-Map ")
    if productVersion or not os.path.exists(circle_bed):
        #Remove MtDNA
        rmMT_cmd = " ".join([bedtools, "intersect -b", config['chr_bed'], "-a", args.sample_bam, ">", sample_bam, "&&", samtools, "index", sample_bam])
        result_note_list.append(rmMT_cmd)
        myc.run_command(rmMT_cmd, ">>>Remove Mt eccDNA Failed")
        myc.rm_files(output_path + '/qname.{}.sorted.bam.tmp*'.format(prefix))
        sort1_cmd = " ".join([samtools, "sort -n -@", cpu_circle, "-o", qname_bam, sample_bam])
        result_note_list.append(sort1_cmd)
        myc.run_command(sort1_cmd, ">>>STAR samtools sort by read name Failed")

        ReadExtractor_cmd = " ".join([config["Circle-Map"], "ReadExtractor -i ", qname_bam, "-o", circular_bam])
        result_note_list.append(ReadExtractor_cmd)
        myc.run_command(ReadExtractor_cmd, ">>>STAR ReadExtractor Failed")
        
        sort2_cmd = " ".join([samtools, "sort -@", cpu_circle, "-o", circular_sort_bam, circular_bam, \
                "&& ", samtools, "index", circular_sort_bam])
        result_note_list.append(sort2_cmd)
        myc.run_command(sort2_cmd, ">>>STAR samtools sort circular bam Failed")
        	
    # readme
    result_note_list.append("End time %s"%time.strftime("%Y-%m-%d %X", time.localtime()))
    myc.write_readme(output_path + "/Pre_Circle-Map.README", version, result_note_list)

    # remove files
    if productVersion:
        myc.rm_files(circular_bam)

