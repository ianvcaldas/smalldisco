#!/bin/bash

## Test script for smalldisco vignette

python ../smalldisco.py sirna -o sirna.bed -a sample.gtf BAM/

python ../smalldisco.py tail -o tails.tsv -g genome.fa --tailor_command ../Tailor/bin/tailor_v1.1_linux_static sirna.bed BAM/
