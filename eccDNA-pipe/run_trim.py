#!/usr/bin/env python

version = "v0.1.0.20220610"

import os,os.path,sys,re
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import myc
import time

if __name__ == '__main__':
    #Parameters to be input.
    parser=ArgumentParser(description="RNA Fastq QC module, version {}".format(version))
    parser.add_argument("-fqs","--fqs", action="store", dest="fqs", 
            help="fastq, for paired end, split by comma", required=True)
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

    # software and args
    trim_galore = config["trim_galore"]
    cpu_trim = config["cpu_trim"]
    mem_trim = int(config["mem_trim"])/1000
    mem_trim = str(int(mem_trim)) + "G"

    # result note
    result_note_list = []
    fqs = args.fqs
    fq_list = fqs.split(",")
    myc.check_files(fq_list)

    # trim_galore Running
    clean_read1 = "{}_1_val_1.fq.gz".format(dirprefix)
    clean_read2 = "{}_2_val_2.fq.gz".format(dirprefix)

    result_note_list.append("Start time %s"%time.strftime("%Y-%m-%d %X", time.localtime()))
    result_note_list.append("Running trim_galore")
    if productVersion or not os.path.exists(clean_read1):
        trim_cmd= " ".join([config["trim_galore"], config["trim_params"], "-o ", output_path, fqs.replace(",", " ")])
        result_note_list.append(trim_cmd)
        myc.run_command(trim_cmd, ">>>Fastq QC Failed")
        	


    # readme
    result_note_list.append("End time %s"%time.strftime("%Y-%m-%d %X", time.localtime()))
    myc.write_readme(output_path + "/trim.README", version, result_note_list)

    # remove files
    if productVersion:
        myc.rm_files("{}*Log*".format(dirprefix))
