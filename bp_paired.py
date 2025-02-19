#!/usr/bin/env python

'''
staatistic paired length
'''

import argparse
import time


def main():
    start_time = time.time()
    
    parser = argparse.ArgumentParser(description='Statictics of paired stem fraction')
    parser.add_argument('-c', required=True, help='Path to the ct file')
    parser.add_argument('-o', help='output file path')
    args = parser.parse_args()

    ct_path = args.c
    out_path = args.o

    with open(ct_path,"r") as f:
        lines = f.readlines()

    seq_start_indices = []
    for i, line in enumerate(lines):
        if "ENST" in line or "chr" in line:
            seq_start_indices.append(i)


    for bp in [6,7,8,9,10]:
        mylog = open(out_path, mode = 'a',encoding='utf-8')
        for seq_idx, seq_start in enumerate(seq_start_indices): 
            if seq_idx == len(seq_start_indices) - 1:
                seq_lines = lines[seq_start:]
            else:
                seq_lines = lines[seq_start:seq_start_indices[seq_idx+1]]
            seq_length = int(seq_lines[0].strip().split()[0])
            bp_pairs = [tuple(map(int, line.split()[4:6])) for line in seq_lines[1:] if int(line.split()[0]) <= seq_length]
            count=1;paired_bp=0;stem_length=[]
            for i in range(len(bp_pairs)):
                if i<(len(bp_pairs)-1) and bp_pairs[i][0]==bp_pairs[i+1][0]+1 and bp_pairs[i][0]>1:
                #让i在倒数第二个停下来，确保是相连的，限制1，0认为是配对的情况
                    count += 1
                else:
                    if count ==bp:
                        paired_bp+=count
                        stem_length.append(count)
                    count=1
            paired_bp_percent = paired_bp / seq_length 
            print(f"{bp}\t{paired_bp}\t{seq_length}\t{paired_bp_percent}\t{stem_length}",file=mylog)
    mylog.close()
        
    end_time = time.time()
    total_time = end_time - start_time
    print("Script finished in {:.2f} seconds".format(total_time))

if __name__ == '__main__':
    main()
