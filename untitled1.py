#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 10:10:41 2024

@author: yamkelamgwatyu
"""

from Bio import SeqIO

def trim_fastq(input_file, output_file, trim_length):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fastq"):
            trimmed_record = record[trim_length:]
            SeqIO.write(trimmed_record, outfile, "fastq")

input_file = "/Users/yamkelamgwatyu/Documents/CSS2024_training/H-clipped_reads.fastq"  # replace with your input FASTQ file path
output_file = "output_trimmed.fastq"  # replace with your desired output file path
trim_length = 500

trim_fastq(input_file, output_file, trim_length)



