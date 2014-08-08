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

        print ">%s/%d_%d" % ( seq_name, s, e) 
        print sseq

