#!/usr/bin/env python
'''
input sequence fasta file and output all k-mers with appindex
'''

import argparse
from collections import Counter
import pandas as pd
import concurrent.futures
import time

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(description='Calculate app_index according to slide windows from sequences')
    parser.add_argument('-s', required=True, help='Path to the FASTA file')
    parser.add_argument('-w', type=int, required=True, help='Size of slide window')
    parser.add_argument('-o', help='output file path')
    args = parser.parse_args()

    fasta_path = args.s
    wsize = args.w
    out_path = args.o

    with open(fasta_path, 'r') as sequences:
        lines = sequences.readlines()

    id_seq={}
    for line in lines:
        if line.startswith('>'):
            genename = line.strip()
            id_seq[genename] = ''
        else:
            id_seq[genename] += line.strip()
    print("fasta file loaded")

    def build_window(sequence,wsize):
        window_seq={}
        n_seq= len(sequence)- wsize +1
        for i in range(n_seq):
            window_seq[i]=sequence[i :i + wsize]
        return window_seq

    mylog = open(out_path, mode = 'a',encoding='utf-8')
    for genename,seq in id_seq.items():
        window_seq=build_window(seq,wsize)
        for key,value in window_seq.items():
            seq = value.replace('C', 'Y').replace('T', 'Y').replace('A', 'R').replace('G', 'R')
            appindex = sum([seq[i] != seq[i+1] for i in range(len(seq)-1)])/(len(seq)-1)
            print(f"{genename}\t{key}\t{appindex}\t{seq}\t{value}",file=mylog)

    mylog.close()
    print("Finished!")

if __name__ == '__main__':
     main()
