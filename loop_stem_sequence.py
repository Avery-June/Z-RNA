#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

parser = argparse.ArgumentParser(description='get the loop or stem sequence/lines from ct file')
parser.add_argument('-c',required=True,help='path to ct files')
parser.add_argument('-b',required=True,type=int,help='min stem length')
parser.add_argument('-o',required=True,help='path to output file folder')
args= parser.parse_args()

ct_path=args.c
stem_length=args.b
out_path=args.o

with open(ct_path,"r") as f:
        lines = f.readlines()

seq_start_indices = []
for i, line in enumerate(lines):
    if "|" in line:
        seq_start_indices.append(i)

sequences = [];loop_sequences=[];stem_lines=[];loop_lines=[]
for seq_idx, seq_start in enumerate(seq_start_indices):
    if seq_idx == len(seq_start_indices) - 1:
        seq_lines = lines[seq_start:]
    else:
        seq_lines = lines[seq_start:seq_start_indices[seq_idx+1]]

    seq_length = int(seq_lines[0].strip().split()[0])
    bp_pairs = [tuple(map(int, line.split()[4:6])) for line in seq_lines[1:] if int(line.split()[5]) <= seq_length]

    count = 1; paired_bp=0; paired_start=[1];seq=[]
    for i in range(len(bp_pairs)):
        if i<(len(bp_pairs)-1) and bp_pairs[i][0]==bp_pairs[i+1][0]+1 and bp_pairs[i][0]>1:
            count += 1
        else:
            if count >=stem_length:
                paired_start.append(i-count+2)
                paired_start.append(i+2)
                seq=seq_lines[i+2-count:i+2]
                stem_lines.extend(seq)
                sequences.append(''.join([line.split('\t')[1] for line in seq]))
            count=1
    paired_start.append(len(seq_lines))


    for idx, start in enumerate(paired_start):
        if  idx % 2==0 and idx<len(paired_start)-1:
            unpaired_lines = seq_lines[paired_start[idx]:paired_start[idx+1]]
            loop_lines.extend(unpaired_lines)
            loop_sequences.append(''.join([line.split('\t')[1] for line in unpaired_lines]))

#The output line is used to find the edit site on the loop
loopline_file = os.path.join(out_path, "loopline.txt")
stemline_file = os.path.join(out_path, "stemline.txt")
loop_seq_file = os.path.join(out_path, "loop_seq.txt")
stem_seq_file = os.path.join(out_path, "stem_seq.txt")

f = open(loopline_file,'w')
f.writelines(loop_lines)
f.close()

f = open(stemline_file,'w')
f.writelines(stem_lines)
f.close()

#output sequence
f = open(loop_seq_file,'w')
for line in loop_sequences:
    f.write(line+'\n')
f.close()

f = open(stem_seq_file,'w')
for line in sequences:
    f.write(line+'\n')
f.close()
