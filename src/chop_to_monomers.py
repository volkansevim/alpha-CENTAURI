#!/usr/bin/env python
from pbcore.io import FastaReader
import re
import sys
import os

hmm_model = sys.argv[1]
in_seq_file = sys.argv[2]

#TODO: use tmp fle than test.tbl and test.out

os.system("rm test_ctrl.tbl; rm test_ctrl.out; nhmmscan --cpu 32 --tblout test_ctrl.tbl -o test_ctrl.out --notextw %s %s" % (hmm_model, in_seq_file)) 


rc_map = dict( zip("ACGTacgtN","TGCAtgcaN") )
seq_db = {}
for r in FastaReader(in_seq_file):
    seq = r.sequence
    seq_db[r.name] = seq
    
hmm_db = {}
with open("test_ctrl.tbl") as f:
    for l in f:
        l = l.strip().split()
        if l[0][0] == "#":
            continue
        seq_name = l[2]
        seq = seq_db[seq_name]
        #print l
        ms = int(l[4])
        me = int(l[5])
        if me - ms < 160:
            continue
        s = int(l[8])
        e = int(l[9])
        #print s, e, e-s
        rev = False
        if s > e:
            #seq = "".join([rc_map[c] for c in seq[::-1]])
            s, e = e, s
            rev = True
        if rev:
            sseq = "".join([rc_map[c] for c in seq[e:s:-1]])
        else:
            sseq = seq[s:e]
        hmm_db.setdefault(seq_name, []).append((s, e))        
        print ">%s/%d_%d" % ( seq_name, s, e) 
        print sseq

flankingout = open ("flanking_seq.fa", 'w')
internalout = open ("interval_seq.fa", 'w')

for seq_name, ranges in hmm_db.items():
    ranges.sort(key = lambda x: x[0])
    beg_seq = 0
    end_seq = len(seq_db[seq_name])
    s_flag = True
    prev_start = 0
    prev_end = end_seq+1
    internalranges = []
    flankingranges = []
    for idx in range(len(ranges)):

        r = ranges[idx]
        if r[0] < prev_end:
            if r[1] > prev_end:
                prev_end = r[1]
        else:
            if s_flag:
                flankingranges.append((prev_end, r[0]))
            else:
                internalranges.append((prev_end, r[0]))
                prev_start = r[0]
                prev_end = r[1]
        s_flag = False
    flankingranges.append((prev_end, end_seq))

    sseq = seq_db[seq_name]
    print >>sys.stderr, internalranges
    print >>sys.stderr, flankingranges
    for irange in internalranges:
        i_s, i_e = irange[0], irange[1] 
        internalout.write(">%s_%i_%i\n%s\n" %(seq_name, i_s, i_e, sseq[i_s:i_e]))
    for frange in flankingranges:     
        f_s, f_e = frange[0], frange[1]
        flankingout.write(">%s_%i_%i\n%s\n" %(seq_name, f_s, f_e, sseq[f_s:f_e]))




