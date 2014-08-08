from falcon_kit import kup, falcon, DWA, get_consensus, get_alignment
import sys

if len(sys.argv) != 3:
    print "Need two arguments: preads_filename inferred_monomer_filename"
    exit(-1)

pread_filename = sys.argv[1]
inferred_monomer_filename = sys.argv[2] 
   
monomer_db = {}
from pbcore.io import FastaReader
seq_db = {}
for r in FastaReader(pread_filename):
#for r in FastaReader("pread_CHM10x.fa"):
    if len(r.sequence) < 8000:
        continue
    seq_db[r.name] = r.sequence
for r in FastaReader(inferred_monomer_filename):
#for r in FastaReader("pread_CHM10x_monomers.fa"):
    #r.name = 87b8fe32_37569_1815/2215_2384
    rid, rng = r.name.split("/")
    if rid not in seq_db:
        continue
    rng = rng.split("_")
    rng = int(rng[0]), int(rng[1])
    monomer_db.setdefault(rid, [])
    monomer_db[rid].append( (rng, r.sequence) )

rc_map = dict(zip("ACGTacgtNn","TGCAtgcaNn"))

import networkx as nx
import numpy as np

for rid, monomers in monomer_db.items():


    aln_data = []
    for i in range(len(monomers)):
        for j in range( i + 1, len(monomers) ):
            t_seq = monomers[i][1]
            q_seq = monomers[j][1]
            alignment = DWA.align(q_seq, len(q_seq), t_seq, len(t_seq), 50, 1)
            aln_str1 = alignment[0].q_aln_str
            aln_str0 = alignment[0].t_aln_str
            aln_size = alignment[0].aln_str_size
            aln_dist = alignment[0].dist
            aln_q_s = alignment[0].aln_q_s
            aln_q_e = alignment[0].aln_q_e
            aln_t_s = alignment[0].aln_t_s
            aln_t_e = alignment[0].aln_t_e
            if aln_size != 0:
                aln_data.append( (i, j, 1.0*aln_dist/aln_size) )
            DWA.free_alignment(alignment)

    for threshold in (0.08, 0.1, 0.12, 0.14):
        G = nx.Graph()
        idt_in_clusters = []
        idt_out_clusters = []
        for i, j, idt in aln_data:
            if idt < threshold:
                G.add_edge(i, j)
                idt_in_clusters.append(1.0-idt)
            else:
                idt_out_clusters.append(1.0-idt)


        c_idx = 0
        data = []
        data_c = {}
        for C in nx.connected_components(G):
            for idx in C:
                s,e = monomers[idx][0]
                data.append( (s, c_idx) )
                data_c.setdefault( c_idx, [])
                data_c[c_idx].append( (s, c_idx) )
            c_idx += 1
        
        if c_idx <= 1: 
            print rid, "NA"
            continue
        
        data.sort()
        new_idx = 0
        idx_map = {}
        for x, y in data:
            if y not in idx_map:
                idx_map[y] = new_idx
            new_idx += 1
   
        x, y = zip(*data)
        y = [ idx_map[c] for c in y ]
        x = np.array(x) 
        l_seq = len(seq_db[rid])
        interval = x[1:] - x[:-1]
        c_intervals = []
        #calcule interval within monomer in a cluster
        for c_index in data_c:
            x = np.array( [ c[0] for c in data_c[c_index] ] )
            x.sort()
            c_interval = np.median( x[1:] - x[:-1] )
            c_intervals.append(c_interval)



        print rid, l_seq, threshold, len(monomers),  len(monomers)*170.0/l_seq, c_idx,\
              round(np.median(c_intervals)/170.0), round(min(c_intervals)/170.0), round(max(c_intervals)/170.0),\
              np.median([round(c) for c in interval/2.0])*2, np.min(interval), np.max(interval),\
              np.mean(idt_in_clusters), np.mean(idt_out_clusters)
