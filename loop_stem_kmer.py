#!/usr/bin/env python
'''
Performs k-mer enrichment analysis by comparing two FASTA datasets (a / b)
'''

import argparse
from collections import Counter
import pandas as pd
import concurrent.futures
import time
import numpy as np

def load_fa(path):
    '''
    a function to load fasta file
    '''
    id_seq = {}
    with open(path, 'r') as sequences:
        lines = sequences.readlines()
    for i,line in enumerate(lines):
        id_seq[i] = line.strip()
    return id_seq

def build_kmerdict(idseq,ksize):
    '''build kmer and statistic numbers '''
    kmers = []#clear the dict
    kmer_freq = []
    for value in idseq.values():       
        n_kmers= len(value)- ksize +1
        
        for i in range(n_kmers):
            kmer=value[i :i + ksize]
            kmers.append(kmer)
        
    kmer_freq=dict(Counter(kmers))
    return kmer_freq

def ratio(freq):
    total_value = sum(freq.values())
    proportion={key: (100 * value / total_value) for key, value in freq.items()}
    return proportion

def fc(dic_ar, dic_br):
    keys = np.intersect1d(list(dic_ar.keys()), list(dic_br.keys()))
    result = {}
    dic_a=ratio(dic_ar)
    dic_b=ratio(dic_br)
    # Calculate ratios for keys that exist in both dictionaries
    for k in keys:
        result[k] = '{:.3f}'.format(dic_a[k] / dic_b[k])
    # Set default value to 1 for keys in dic_a that are not in dic_b
    for k in set(dic_ar.keys()) - set(dic_br.keys()):
        result[k] = '{:.3f}'.format(dic_a[k] / (1.0/sum(dic_br.values())))

    # Set default value to 1 for keys in dic_b that are not in dic_a
    for k in set(dic_b.keys()) - set(dic_a.keys()):
        result[k] = '{:.3f}'.format( (1.0/sum(dic_ar.values())) / dic_b[k])

    return result


def app_index(dict):
    result=[]
    for key,value in dict.items():
        seq = key.replace('C', 'Y').replace('U', 'Y').replace('T', 'Y').replace('A', 'R').replace('G', 'R')
        appindex ='{:.3f}'.format(sum([seq[i] != seq[i+1] for i in range(len(seq)-1)])/(len(seq)-1))
        output = (appindex,value,key,seq)
        result.append(output)
    return result

def main():
    start_time = time.time()
    
    parser = argparse.ArgumentParser(description='Process FASTA file and k-mer size')
    parser.add_argument('-a', required=True, help='Path to the FASTA file')
    parser.add_argument('-b', required=True, help='Path to the FASTA file')
    parser.add_argument('-ksize', type=int, required=True, help='Size of k-mer')
    parser.add_argument('-o', help='output file path')
    args = parser.parse_args()

    afasta_path = args.a
    bfasta_path = args.b
    ksize = args.ksize
    path = args.o

    load_start_time = time.time()
    aid_seq = load_fa(afasta_path)
    print("Fasta A loaded")
    bid_seq = load_fa(bfasta_path)
    print("Fasta B loaded")
    load_end_time = time.time() 
    print("Fasta files loaded in {:.2f} seconds".format(load_end_time - load_start_time))

    print("Building k-mer dictionaries...")
    build_start_time = time.time()
    with concurrent.futures.ThreadPoolExecutor() as executor:
        akmer_future = executor.submit(build_kmerdict, aid_seq, ksize)
        bkmer_future = executor.submit(build_kmerdict, bid_seq, ksize)
        
        akmer_freq = akmer_future.result()
        print("A k-mer dictionary built")
        
        bkmer_freq = bkmer_future.result()
        print("B k-mer dictionary built")
    build_end_time = time.time()
    print("kmer dict built in {:.2f} seconds".format(build_end_time - build_start_time)) 
    
    print("Calculating fold changes...")
    #start_time = time.time()
    foldchange = fc(akmer_freq, bkmer_freq)
    #print("Fold changes calculated")
    end_time = time.time()
    
    print("Calculating app indexes...")
    #start_time = time.time()
    app_index_2 = app_index(foldchange)
    
    print("Saving results to file...")
    app_index_df = pd.DataFrame(app_index_2)
    app_index_df.columns = ["app_index", "foldchange",  "kmer", "kmer(YR)"]
    app_index_df.to_csv(path, sep='\t', index=False, header=True)
    print("Results saved to", path)
    
    end_time = time.time()
    total_time = end_time - start_time
    print("Script finished in {:.2f} seconds".format(total_time))

if __name__ == '__main__':
    main()
