#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 14:00:38 2020

@author: janbinkowski
"""
import os
import subprocess

def alignment(reference: str, read_1: str, read_2: str, raw_project_dir: str) -> str:
    
    alignment_report_path = os.path.join(raw_project_dir, "aligment_raport.txt")
    sam_file_path = os.path.join(raw_project_dir, "raw_alignment.sam")
    
    querry = ["hisat2", "-x", reference, "-1", read_1, "-2", read_2, 
              "--summary-file", alignment_report_path, 
              "--no-spliced-alignment", "-S", sam_file_path]

    subprocess.run(querry, text=True)
    return sam_file_path

def sam_to_bam(raw_project_dir: str, sam_file_path: str) -> None:
    
    bam_file = os.path.join(raw_project_dir, "raw_alignment.bam")
    querry = ["samtools", "view", sam_file_path, "-b", "-o", bam_file]

    subprocess.run(querry)
    
def slop(regions_bed: str, genome_ref: str, output: str) -> str:
    
    outputfile = os.path.join(output, "slop_500bp.bed")
    querry = ["bedtools", "slop", "-i", regions_bed, "-g", genome_ref, "-b", 500, ">", outputfile]
    subprocess.run(querry)
    
    return outputfile


    