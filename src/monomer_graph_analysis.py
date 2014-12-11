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

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('reads_fasta',
                    action='store',
                    help='Name of the FASTA file containing the reads.')

parser.add_argument('monomer_fasta',
                    action='store',
                    help='Name of the inferred monomer file.')

parser.add_argument('-l',
                    action='store',
                    type=int,
                    default=170,
                    dest='monomer_length',
                    help='Average length of a monomer.')

parser.add_argument('-d',
                    action='store',
                    type=int,
                    default=5,
                    dest='head_to_tail',
                    help='Maximum allowed head-to-tail distance \
                          between two adjacent monomers.')

parser.add_argument('-s',
                    action='store',
                    type=int,
                    default=2000,
                    dest='shortest_read_len',
                    help='Minimum allowed read length.')

parser.add_argument('-t',
                    action='append',
                    type=float,
                    dest='thresholds_list',
                    default=DefaultList([
                        0.98, 0.97, 0.96, 0.95, 0.94, 0.93,
                        0.92, 0.91, 0.9, 0.89, 0.88
                        ]),
                    help='Specifies a clustering threshold. Multiple -t\
                    allowed. Values sorted and tested in descending\
                    order.',
                    )

parser.add_argument('--version', action='version', version='%(prog)s 0.2')

results = parser.parse_args()
pread_filename = results.reads_fasta
inferred_monomer_filename = results.monomer_fasta
len_threshold = results.shortest_read_len
average_monomer_len = results.monomer_length
identity_thresholds = results.thresholds_list
allowed_max_head_to_tail = results.head_to_tail

header_items = ("RID",
                "Regularity",
                "Read_Len",
                "Thresh",
                "#All_Monomers(clustered+not_clustered)",
                "#mono_in_a_cluster",
                "Isolates_(unclustered_monomers)",
                "Clustered_monomer_fraction_in_read",
                "#total_clusters_(distinct_monomers_in_HOR)",
                "mean_identity_within_clusters",
                "mean_identity_between_clusters",
                "min_monomeric_period",
                "max_monomeric_period",
                "median_monomeric_period",
                "Normalized_min_monomeric_period",
                "Normalized_max_monomeric_period",
                "min_head_to_tail_interval",
                "max_head_to_tail_interval",
                "median_head_to_tail_interval"
                )
header = "\t".join(header_items)
monomer_db = {}
seq_db = {}
rc_map = dict(zip("ACGTacgtNn", "TGCAtgcaNn"))
stats_formatting = "".join((
    "%d\t%f\t%d\t%d\t%d\t%f\t%d\t%f\t%f\t",
    "%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\n"))
filename_root = os.path.basename(pread_filename).split(".")[0]
regular_HORs_file = open(filename_root + "_regularHORs.fa", 'w')
inversions_file = open(filename_root + "_inversions.fa", 'w')
irregular_HORs_file = open(filename_root + "_irregularHORs.fa", 'w')
too_short_reads_file = open(filename_root + "_too_short_reads.fa", 'w')
no_HOR_reads_file = open(filename_root + "_no_HOR_reads.fa", 'w')
regular_pattern_file = open(filename_root + "_regularHORs_pattern.txt", 'w')
irregular_pattern_file = \
    open(filename_root + "_irregularHORs_pattern.txt", 'w')
inversions_pattern_file = open(filename_root + "_inversions_pattern.txt", 'w')
stats_file = open(filename_root + "_stats.txt", 'w')
stats_file.write(header + "\n")

# Print parameters
print "Average monomer length: ", average_monomer_len
print "Max head-to-tail distance: ", allowed_max_head_to_tail
print "Shortest read length: ", len_threshold
print "Clustering thresolds: ", identity_thresholds

# IMPORT FASTA FILES #
for r in FastaReader(pread_filename):
    # Load all reads from the pread file seq_db.
    # seq_db[Read_ID] = sequence
    if len(r.sequence) < len_threshold:
        too_short_reads_file.write(">"+r.name + "\n" + r.sequence + "\n")
        continue
    seq_db[r.name] = r.sequence

# Load all monomers found in preads into monomer_db.
# monomer_db[Read_ID] = [(start, end), sequence]
for r in FastaReader(inferred_monomer_filename):
    # Parse the read tag.
    # Tag Format: ReadID/RangeStart_RangeEnd/Orientation
    rid, rng, orientation = r.name.split("/")
    # Skip if the read doesn't have any monomers.
    if rid not in seq_db:
        continue
    rng = rng.split("_")
    rng = int(rng[0]), int(rng[1])
    monomer_db.setdefault(rid, [])
    monomer_db[rid].append((rng, r.sequence, orientation))

print len(seq_db), " sequences read.",\
    "Reads withmonomers:", len(monomer_db.keys())

# RUN OVER ALL READS THAT CONTAIN MONOMERS #
for rid, monomers in monomer_db.items():
    aln_data = []
    range_list = []
    total_monomer_len = 0
    # Align all monomers on the read to each other.
    for i in range(len(monomers)):
        for j in range(i + 1, len(monomers)):
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

            # Skip if no alignment
            if aln_size != 0:
                # Store the alignment score for the monomer pair: (i, j, score)
                aln_data.append((i, j, 1 - (1.0*aln_dist/aln_size)))
            DWA.free_alignment(alignment)

        mono_range = monomers[i][0]
        # Add up app monomer ranges to calculate total monomer content in read
        total_monomer_len += (mono_range[1] - mono_range[0])
        # Store start-end coordinates of the monomer.
        range_list.append(mono_range[0])
        range_list.append(mono_range[1])

    is_regular = False
    inversion_detected = False

    # SCAN THE READ FOR AN HOR AT MULTIPLE CLUSTERING THRESHOLDS #
    for threshold in identity_thresholds:
        G = nx.Graph()
        idt_in_clusters = []
        idt_out_clusters = []
        # Run over each monomer pair in aln_data to cluster monomers.
        for i, j, idt in aln_data:
            # Connect monomers i and j, if score is larger than threshold.
            # Monomers in the same cluster are treated as the same kind.
            if idt >= threshold:
                G.add_edge(i, j)
                idt_in_clusters.append(idt)
            else:
                idt_out_clusters.append(idt)

        cluster_count = 0
        data = []
        data_c = {}
        l_seq = len(seq_db[rid])
        forward = False
        reverse = False

        # Check if monomers in either orientation exist in the read.
        for mono in monomers:
            orientation = mono[2]
            if orientation == "F":
                forward = True
            else:
                reverse = True
        # Mark read as INVERTED if monomers in both orientations exist.
        if forward and reverse:
            inversion_detected = True

        # Run over all clusters
        for C in nx.connected_components(G):
            # Run over all nodes (monomers) in cluster.
            for idx in C:
                s, e = monomers[idx][0]
                # Store (monomer_start, end, cluster_index) in 'data.'
                # Store (monomer_start, cluster_index) in 'data_c.'
                data.append((s, e, cluster_count))
                data_c.setdefault(cluster_count, [])
                data_c[cluster_count].append((s, cluster_count))
            cluster_count += 1

        if cluster_count <= 1 and not inversion_detected:
            # No or only one cluster. No HOR detected.\
            # Test with a different threshold.
            continue

        # Analyze clusters, if more than one detected.
        if cluster_count > 1:
            # Sort detected monomers by their coordinates.
            data.sort()
            # Create a symbolic HOR pattern, i.e., ABCDABCDABCD
            symbolic_pattern = ""
            for s, e, c in data:
                symbolic_pattern += chr(65 + c)

            # Calculate head-to-tail distances between clustered monomers.
            x, e, y = zip(*data)
            x = np.array(x)
            head_to_head_intervals = x[1:] - x[:-1]
            head_to_tail_intervals = x[1:] - e[:-1] - 1

            # Calculate intervals between monomers in the same cluster
            c_intervals = []
            all_monomer_periods = []
            for c_index in data_c:
                x = np.array([c[0] for c in data_c[c_index]])
                x.sort()
                # Periods: intervals between monomers of the same type
                monomer_periods = x[1:] - x[:-1]
                all_monomer_periods.extend(monomer_periods)

            min_monomer_period = min(all_monomer_periods)
            max_monomer_period = max(all_monomer_periods)
            median_monomer_period = round(np.median(all_monomer_periods))
            min_head_to_tail = min(head_to_tail_intervals)
            max_head_to_tail = max(head_to_tail_intervals)
            median_head_to_tail = round(np.median(head_to_tail_intervals))
            max_abs_head_to_tail = max(abs(head_to_tail_intervals))
            isolate_count = len(monomers) - len(data)
            HOR_start = min(range_list)
            HOR_end = max(range_list)
            monomeric_fraction_in_HOR = \
                1.0*total_monomer_len/(HOR_end - HOR_start)

            # Normalize monomer periods by the median.
            if median_monomer_period > 0:
                normalized_min_monomer_period = \
                    min_monomer_period/median_monomer_period
                normalized_max_monomer_period = \
                    max_monomer_period/median_monomer_period
            else:
                # Assign two out-of-range numbers if median is zero.
                normalized_min_monomer_period = 0
                normalized_max_monomer_period = 2

        # Exit the threshold loop if regularity or inversion detected\
        if inversion_detected:
            break
        elif ((max_abs_head_to_tail <= allowed_max_head_to_tail) and
                (isolate_count == 0) and
                (normalized_max_monomer_period <= 1.05) and
                (normalized_min_monomer_period >= 0.95) and
                (not inversion_detected)):
            # Mark as regular
            is_regular = True
            break
    # END OF THE THRESHOLD LOOP #

    # PREPARE THE OUTPUT BUFFER #
    if cluster_count > 1:
        # Create a new fasta tag:
        # >ReadID___ReadLen__HORstart_HORend__HORperiod
        fasta_tag = "".join((
            ">",
            rid,
            "___",
            str(len(seq_db[rid])),
            "__",
            str(HOR_start),
            "_",
            str(HOR_end),
            "__HOR",
            str(int(cluster_count)),
            "\n"))
        fasta_seq = "".join((seq_db[rid][HOR_start:HOR_end+1], "\n"))

        # Create the output buffer for regular class.
        stats = stats_formatting % \
            (l_seq,
             threshold,
             len(monomers),
             len(data),
             isolate_count,
             1.0*len(data)*average_monomer_len/l_seq,
             cluster_count,
             np.mean(idt_in_clusters),
             np.mean(idt_out_clusters),
             min_monomer_period, max_monomer_period,
             median_monomer_period,
             normalized_min_monomer_period,
             normalized_max_monomer_period,
             min_head_to_tail,
             max_head_to_tail,
             median_head_to_tail)
    else:
        # cluster_count <= 1. Inversion detected, but not an HOR
        # Create the output buffer for inversion class.
        fasta_tag = ">" + rid + "\n"
        fasta_seq = seq_db[rid] + "\n"
        stats = "%d %f %d %d %d %s %d %s\n" %\
            (l_seq,
             threshold,
             len(monomers),
             len(data),
             len(monomers),
             ".",
             cluster_count,
             "-1 "*10)
        symbolic_pattern = "N/A (No HOR)"

    # WRITE OUT THE RESULTS #
    if inversion_detected:
        # Read with inversion.
        stats_file.write(rid + " V " + stats)
        inversions_file.write(fasta_tag)
        inversions_file.write(fasta_seq)
        inversions_pattern_file.write(fasta_tag)
        inversions_pattern_file.write(symbolic_pattern + "\n")
    elif is_regular:
        # A read with a regular HOR.
        stats_file.write(rid + " R " + stats)
        regular_HORs_file.write(fasta_tag)
        regular_HORs_file.write(fasta_seq)
        regular_pattern_file.write(fasta_tag)
        regular_pattern_file.write(symbolic_pattern + "\n")
    elif cluster_count > 1:
        # An irregular HOR.
        # (Some monomers clustered, but monomer order is irregular)
        stats_file.write(rid + " I " + stats)
        irregular_HORs_file.write(fasta_tag)
        irregular_HORs_file.write(fasta_seq)
        irregular_pattern_file.write(fasta_tag)
        irregular_pattern_file.write(symbolic_pattern + "\n")
    else:
        # No HOR detected at the minimum identity threshold.
        # (There is no more than one cluster.)
        stats_file.write(rid + " N " + stats)
        no_HOR_reads_file.write(fasta_tag)
        no_HOR_reads_file.write(fasta_seq)
# END OF THE READ LOOP

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
