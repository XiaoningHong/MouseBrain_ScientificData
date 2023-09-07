#!/usr/bin/env python
#-*- coding: utf-8 -*-
# PROGRAM : workflow
# AUTHOR  : xiaoning.hong@outlook.com
# CREATED : Thu June 09 2022
# VERSION : V1.0.0
# UPDATE  : [V1.0.0] Thu June 09 2022

import os,sys,re
from argparse import ArgumentParser
from collections import defaultdict
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '/home/hongxiaoning/pipeline/bin'))
import myc
from pyflow import WorkflowRunner

def sample_list(sample_list):
    sample_dict = {}
    with open(sample_list, 'r') as f_sample:
        for line in f_sample:
            ls  = re.split("\t", line.strip())
            if re.match("#.*", line) or ls[0] == 'sample':
                continue
            sample_dict[ls[0]] = [ls[1], ls[2]]

    return sample_dict

def clean_fq(project_dir, sample, fq1, fq2):
    read1_name = os.path.basename(fq1).replace('.fq','_val_1.fq.gz')
    read1_name = os.path.basename(fq1).replace('.fq.gz','_val_1.fq.gz')
    read1_name = os.path.basename(fq1).replace('.fastq','_val_1.fq.gz')
    read1_name = os.path.basename(fq1).replace('.fastq.gz','_val_1.fq.gz')
                
    read2_name = os.path.basename(fq2).replace('.fq','_val_2.fq.gz')
    read2_name = os.path.basename(fq2).replace('.fq.gz','_val_2.fq.gz')
    read2_name = os.path.basename(fq2).replace('.fastq','_val_2.fq.gz')
    read2_name = os.path.basename(fq2).replace('.fastq.gz','_val_2.fq.gz')
                
    clean_read1 = os.path.join(project_dir, sample, "trim", read1_name)
    clean_read2 = os.path.join(project_dir, sample, "trim", read2_name)
    clean_fqs = '%s,%s'%(clean_read1, clean_read2)

    return(clean_fqs)

class general_process():
    def __init__(self,params):
        self.params=params
        self.debug = ""
        if self.params.debug:
            self.debug = "--debug"
        configs=self.params.config.split(",")
        self.sample_cfg = configs[0]
    def update_dependence(self,ori_list,update_list,update_type='add'):
        tmp_includes=set(ori_list)
        if update_type=='add':
            tmp_includes.update(update_list)
        if update_type == 'remove':
            tmp_includes.remove(update_list)
        return list(tmp_includes)
    def analysis_step(self):
        all_steps=set(['trim', 'aln', 'circle'])
        includes=[];excludes=[]
        if self.params.exrun!=None:
            exrun=self.params.exrun.strip().lower().split(",")
            for each in exrun:
                if re.match("\+.*",each):
                    includes.append(each.replace("+",""))
                else:
                    excludes.append(each)
        if len(includes) != 0:
            all_steps = set(includes)
        run_steps = all_steps - set(excludes)
        return run_steps

class runPipeline(WorkflowRunner, general_process):
    def __init__(self, params):
        #self.args = args
        general_process.__init__(self,params)

    def update_dependence(self,ori_list,update_list,update_type='add'):
        tmp_includes=set(ori_list)
        if update_type=='add':
            tmp_includes.update(update_list)
        if update_type == 'remove':
            tmp_includes.remove(update_list)
        return list(tmp_includes)

    def workflow(self):
        run_steps = self.analysis_step()        
        all_dependence = []
        pre_dependence=[]

        sample_dict = sample_list(self.params.sample_list)
        myc.check_files([self.params.config])
        config = myc.resolveConfig(self.params.config)
        # software or args
        python        = config["python"]

        #createDirectory
        print ("Creating Directory")
        project_dir = os.path.abspath(self.params.output_dir)
        myc.createDirectory(project_dir)
            
        #Fastq QC , low quality remove and adapter trim
        cpu_trim = int(config["cpu_trim"])
        mem_trim = int(config["mem_trim"])
        if 'trim' in run_steps:
            print ("Start Fastq QC")
            trim_dir = os.path.join(project_dir, "00.trim")
            myc.createDirectory(trim_dir)
            for sample in sample_dict:
                sample_Fastqc_dir = os.path.join(trim_dir, sample)
                raw_fqs = '%s,%s'%(sample_dict[sample][0],sample_dict[sample][1])
                trim_app= config['trim_app']
                trim_cmd = " ".join([python, trim_app, "--fqs", raw_fqs, "--cfg", self.params.config,
                    "--output_path", sample_Fastqc_dir, "--prefix", sample, self.debug])
                print (trim_cmd)
                self.addTask("I_%s_trim"%sample, trim_cmd, nCores=cpu_trim, memMb=mem_trim)

        #Alignment by STAR
        if 'aln' in run_steps:
            print ("Start aln fastq to genome")
            aln_dir = os.path.join(project_dir, "01.bwa")
            myc.createDirectory(aln_dir)
            for sample in sample_dict:
                sample_aln_dir = os.path.join(aln_dir, sample)
                cpu_aln = int(config["cpu_aln"]);mem_aln = int(config["mem_aln"])
                sample_bam = "{}/{}.sort.bam".format(sample_aln_dir, sample)
                markdup="--markDuplicate " + config["markDuplicate"]
                aln_app=config['aln_app']
                if 'trim' not in run_steps:
                    clean_fqs = '%s,%s'%(sample_dict[sample][0],sample_dict[sample][1])
                    aln_dependence=[]
                else:
                    clean_fqs = clean_fq(project_dir, sample, sample_dict[sample][0], sample_dict[sample][1])
                    aln_dependence=["I_%s_trim"%sample]
                
                aln_cmd = " ".join([python, aln_app, "--fqs", clean_fqs, "--cfg", self.params.config,"--output_path", sample_aln_dir, "--prefix", sample, markdup, self.debug])
                print (aln_cmd)
                self.addTask("II_%s_Aln"%sample, aln_cmd, nCores = cpu_aln,memMb = mem_aln,dependencies = aln_dependence)
                pre_dependence = self.update_dependence(pre_dependence, ["II_%s_Aln"%sample], 'add')
        #Count gene and transcript expression (RPKM/TPM)
        if 'circle' in run_steps:
            print ("Start Detecting Circle analysis")
            aln_dir = os.path.join(project_dir, "01.bwa")
            circle_dir = os.path.join(project_dir, "02.CircleMap")
            myc.createDirectory(circle_dir)
            for sample in sample_dict:
                sample_circle_dir = os.path.join(circle_dir, sample)
                if 'aln' in run_steps:
                    pre_circle_dependence = ["II_%s_Aln"%sample]
                else:
                    pre_circle_dependence = []
                circle_cpu = int(config["cpu_circle"])
                circle_mem = int(config["mem_circle"])
                
                if config['markDuplicate'] == "no":
                    sample_bam = "{}/{}/{}.sort.bam".format(aln_dir, sample, sample)
                else:
                    sample_bam = "{}/{}/{}.rmdup.bam".format(aln_dir, sample, sample)
                
                Pre_CircleMap_cmd = " ".join([python, config['Circle-Map_pre'], "--sample_bam", sample_bam, "--cfg", args.config, "--output_path", sample_circle_dir, "--prefix", sample ])
                print (Pre_CircleMap_cmd)
                self.addTask("III_%s_CM_pre"%sample, Pre_CircleMap_cmd, nCores = 4, memMb = circle_mem, dependencies = pre_circle_dependence)
                circle_dependence = self.update_dependence(pre_circle_dependence, ["III_%s_CM_pre"%sample], 'add')
                run_CircleMap_cmd = " ".join([python, config['Circle-Map_app'], "--sample_bam", sample_bam, "--cfg", args.config, "--output_path", sample_circle_dir, "--prefix", sample ])
                print (run_CircleMap_cmd)
                self.addTask("IV_%s_CM"%sample, run_CircleMap_cmd, nCores = circle_cpu, memMb = circle_mem, dependencies = circle_dependence)

if __name__ == "__main__":
    parser = ArgumentParser(description='eccDNA Detecting pipeline')
    parser.add_argument("-c", "--config", action="store", dest="config",
                            help="project special configure file, separate by comma", required=True)
    parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                            help="project output dir", required=True)
    parser.add_argument("-s", "--sample_list", action="store", dest="sample_list",
                            help="sample list", required=True)
    parser.add_argument("--exrun", action="store", dest="exrun",
                            help="the steps you do not wanna perform. \
                            options:trim,aln,circle", default=None)  
    parser.add_argument("--debug", action="store_true", dest="debug",
                            help="for debug module, skipping existed file, default False", default=False)
    args = parser.parse_args()

    wflow = runPipeline(args)
    retval = wflow.run(mode="sge", nCores=200, isContinue="Auto", isForceContinue=False, dataDirRoot=args.output_dir, retryMax=2)
    sys.exit(retval)
