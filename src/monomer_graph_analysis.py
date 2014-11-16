#!/usr/bin/env python
from falcon_kit import kup, falcon, DWA, get_consensus, get_alignment
from pbcore.io import FastaReader
import sys
import re
import networkx as nx
import numpy as np
import os.path
import argparse

class DefaultList(list):
    def __copy__(self):
        return [] 

parser = argparse.ArgumentParser()

parser.add_argument('-f', 
                    action='store', 
                    dest ='fasta_file',
                    help ='Name of the FASTA file containing the reads.',
                    required = True)

parser.add_argument('-m', action = 'store', 
                    dest = 'monomers_file',
                    help = 'Name of the inferred monomer file.',
                    required = True)

parser.add_argument('-l', 
                    action = 'store',
                    type = int,
                    default = 170,
                    dest = 'monomer_length',
                    help = 'Average length of a monomer.')

parser.add_argument('-d', 
                    action = 'store',
                    type = int,
                    default = 5,
                    dest = 'head_to_tail',
                    help = 'Maximum allowed head-to-tail distance between two adjacent monomers.')

parser.add_argument('-s', 
                    action = 'store',
                    type = int,
                    default = 2000,
                    dest = 'shortest_read_len',
                    help = 'Minimum allowed read length.')

parser.add_argument('-t', 
                    action='append', 
                    type=int,
                    dest='thresholds_list',
                    default = DefaultList([0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9, 0.89, 0.88]),
                    help='Specifies a clustering threshold. Multiple -t allowed. Values sorted and tested in descending order.',
                    )

parser.add_argument('--version', action='version', version='%(prog)s 0.2')

results = parser.parse_args()
pread_filename = results.fasta_file
inferred_monomer_filename = results.monomers_file 
len_threshold = results.shortest_read_len
average_monomer_len = results.monomer_length
identity_thresholds = results.thresholds_list
allowed_max_head_to_tail = results.head_to_tail

header = "RID\t" + "Regularity\t" + "Read_Len\t" + "Thresh\t" + "#All_Monomers(clustered+not_clustered)\t" +\
    "#mono_in_a_cluster\t" + "Isolates_(all_minus_clustered)\t" + \
    "Clustered_monomer_fraction_in_read_(HOR%_in_read)\t" + \
    "#total_clusters_(distinct_monomers_in_HORs)\t" + "median_#monomers_in_HOR\t" + \
    "min_#monomers_in_HOR\t" + "max_#monomers_in_HOR\t" + "median_monomer_len\t" + "min_monomer_len\t" + \
    "max_monomer_len\t" + "mean_identity_within_clusters\t" + "mean_identity_between_clusters\t" +\
    "min_monomeric_period\t" + "max_monomeric_period\t" + "median_monomeric_period\t" +\
    "Normalized_min_monomeric_period\t" + "Normalized_max_monomeric_period\t" +\
    "min_head_to_tail_interval\t" + "max_head_to_tail_interval\t" + "median_head_to_tail_interval"
   
monomer_db = {}
seq_db = {}

filename_root = os.path.basename(pread_filename).split(".")[0]
regular_HORs_file = open (filename_root + "_regularHORs.fa", 'w')
inversions_file = open (filename_root + "_inversions.fa", 'w')
irregular_HORs_file = open (filename_root + "_irregularHORs.fa", 'w')
too_short_reads_file = open (filename_root + "_too_short_reads.fa", 'w')
no_HOR_reads_file = open (filename_root + "_no_HOR_reads.fa", 'w')
regular_pattern_file = open (filename_root + "_regularHORs_pattern.txt", 'w')
irregular_pattern_file = open (filename_root + "_irregularHORs_pattern.txt", 'w')
inversions_pattern_file = open (filename_root + "_inversions_pattern.txt", 'w')
stats_file = open (filename_root + "_stats.txt", 'w')
stats_file.write(header + "\n")

for r in FastaReader(pread_filename):
    # Load all reads from the pread file seq_db.
    # seq_db[Read_ID] = sequence
    if len(r.sequence) < len_threshold:
        too_short_reads_file.write(">"+r.name + "\n" + r.sequence + "\n")        
        continue
    seq_db[r.name] = r.sequence

for r in FastaReader(inferred_monomer_filename):
    # Load all monomers found in preads into monomer_db.
    # monomer_db[Read_ID] = [(start, end), sequence]
    rid, rng, direction = r.name.split("/")
    if rid not in seq_db:
        continue
    rng = rng.split("_")
    rng = int(rng[0]), int(rng[1])
    monomer_db.setdefault(rid, [])
    monomer_db[rid].append( (rng, r.sequence, direction) )

#for key in monomer_db.keys():
    #print key
#sys.exit()

print len(seq_db), " sequences read.", "Reads with monomers: ", len(monomer_db.keys())
rc_map = dict(zip("ACGTacgtNn","TGCAtgcaNn"))

for rid, monomers in monomer_db.items():    
    aln_data = []    
    range_list = []
    total_monomer_len = 0
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
                aln_data.append( (i, j, 1 - (1.0*aln_dist/aln_size)) )
            DWA.free_alignment(alignment)            

        mono_range = monomers[i][0]
        total_monomer_len += (mono_range[1] - mono_range[0])
        range_list.append(mono_range[0])
        range_list.append(mono_range[1])
          
    is_regular = False 
    inversion_detected = False
    # Scan the read for an HOR at multiple clustering thresholds 
    for threshold in identity_thresholds:  
        G = nx.Graph()
        idt_in_clusters = []
        idt_out_clusters = []
        for i, j, idt in aln_data:
            if idt >= threshold:
                G.add_edge(i, j)
                idt_in_clusters.append(idt)
            else:            
                idt_out_clusters.append(idt)

        c_idx = 0
        data = []
        data_c = {}
        l_seq = len(seq_db[rid])
        forward = False
        reverse = False
        for C in nx.connected_components(G):
            forward = False
            reverse = False
            for idx in C:
                s,e = monomers[idx][0]
                direction = monomers[idx][2]
                data.append( (s, e, c_idx) )
                data_c.setdefault( c_idx, [])
                data_c[c_idx].append( (s, c_idx) )
                if direction == "F":
                    forward = True
                else:
                    reverse = True                
            c_idx += 1
            if forward and reverse and not inversion_detected:
                inversion_detected = True
        
        if c_idx <= 1 and not inversion_detected: 
            # 0 or only 1 cluster. No HOR detected.     
            continue
        
        data.sort()
        # Create symbolic HOR pattern
        symbolic_pattern = ""
        for s, e, c in data:
            symbolic_pattern += chr(65 + c)
        
        # Calculate monomeric lenghts (variable interval)
        new_idx = 0
        idx_map = {}
        for x, e, y in data:
            if y not in idx_map:
                idx_map[y] = new_idx
            new_idx += 1
        
        x, e, y = zip(*data)
        y = [ idx_map[c] for c in y ]
        x = np.array(x)     
        interval = x[1:] - x[:-1]
        #Calculate Head to Tail distances between clustered monomers
        head_to_tail_intervals =  x[1:] - e[:-1] - 1
        
        # Calculate intervals between monomers in the same cluster
        c_intervals = []    
        all_monomer_periods = []
        for c_index in data_c:
            x = np.array( [ c[0] for c in data_c[c_index] ] )
            x.sort()
            monomer_periods = x[1:] - x[:-1]
            # c_intervals contains the periods of each monomer within a cluster
            c_interval = np.median(monomer_periods)
            c_intervals.append(c_interval)
            all_monomer_periods.extend(monomer_periods)
        
        min_monomer_period = min(all_monomer_periods) 
        max_monomer_period = max(all_monomer_periods) 
        median_monomer_period = round(np.median(all_monomer_periods))            
        min_head_to_tail =  min(head_to_tail_intervals)
        max_head_to_tail =  max(head_to_tail_intervals)
        median_head_to_tail = round(np.median(head_to_tail_intervals))        
        max_abs_head_to_tail = max(abs(head_to_tail_intervals))
        isolate_count = len(monomers) - len(data)
        min_monomers_in_HOR = round(min(c_intervals)/average_monomer_len)
        max_monomers_in_HOR = round(max(c_intervals)/average_monomer_len)
        median_monomers_in_HOR = round(np.median(c_intervals)/average_monomer_len)
        HOR_start = min(range_list)
        HOR_end = max(range_list) 
        monomeric_fraction_in_HOR = 1.0*total_monomer_len/(HOR_end - HOR_start)    
        
        if median_monomer_period > 0:
            normalized_min_monomer_period = min_monomer_period/median_monomer_period
            normalized_max_monomer_period = max_monomer_period/median_monomer_period
        else:
            normalized_min_monomer_period = 0
            normalized_max_monomer_period = 2                       
                
        fasta_tag = ">" + rid + "___" + str(len(seq_db[rid])) + "__" + \
                    str(HOR_start) + "_" + str(HOR_end) + "__HOR" + str(int(min_monomers_in_HOR))+"\n"
        fasta_seq = seq_db[rid][HOR_start:HOR_end+1] + "\n"

        stats =  "%d\t%f\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\n" % \
                    (l_seq, threshold, len(monomers), len(data), isolate_count, len(data)*average_monomer_len/l_seq, c_idx,\
                    median_monomers_in_HOR, min_monomers_in_HOR, max_monomers_in_HOR,\
                    np.median([round(c) for c in interval/2.0])*2, np.min(interval), np.max(interval),\
                    np.mean(idt_in_clusters), np.mean(idt_out_clusters), \
                    min_monomer_period, max_monomer_period, median_monomer_period, \
                    normalized_min_monomer_period, normalized_max_monomer_period, \
                    min_head_to_tail, max_head_to_tail, median_head_to_tail )        
        
        if inversion_detected:
            stats_file.write(rid + " V " + stats)
            inversions_file.write(fasta_tag)
            inversions_file.write(fasta_seq)
            inversions_pattern_file.write(fasta_tag)
            inversions_pattern_file.write(symbolic_pattern + "\n")
            break 
        elif (max_monomers_in_HOR == min_monomers_in_HOR) and \
            (max_abs_head_to_tail <= allowed_max_head_to_tail) and \
            (isolate_count == 0) and (normalized_max_monomer_period <= 1.05) and \
            (normalized_min_monomer_period >= 0.95) and (not inversion_detected): 
            is_regular = True
            # Write out the stats
            stats_file.write(rid + " R " + stats)
            regular_HORs_file.write(fasta_tag)
            regular_HORs_file.write(fasta_seq)            
            regular_pattern_file.write(fasta_tag)
            regular_pattern_file.write(symbolic_pattern + "\n")            
            break
    
    # HOR is irregular
    if not is_regular and not inversion_detected:
        if c_idx > 1 :
            # An irregular HOR detected
            stats_file.write(rid + " I " + stats)
            irregular_HORs_file.write(fasta_tag)
            irregular_HORs_file.write(fasta_seq)
            irregular_pattern_file.write(fasta_tag)
            irregular_pattern_file.write(symbolic_pattern + "\n")
        else: 
            # No HOR detected at the minimum identity threshold
            buf = "%d %f %d %d %d %s %d %s\n" %\
                (l_seq, threshold, len(monomers), len(data), len(monomers), ".", c_idx, "-1 "*16)
            stats_file.write(rid + " N " + buf)
            no_HOR_reads_file.write(">" + rid + "\n" + seq_db[rid] + "\n")  

regular_HORs_file.close()
irregular_HORs_file.close()
inversions_file.close()
no_HOR_reads_file.close()
too_short_reads_file.close()
regular_HORs_file.close()
regular_pattern_file.close()
irregular_pattern_file.close()
inversions_pattern_file.close()
stats_file.close()