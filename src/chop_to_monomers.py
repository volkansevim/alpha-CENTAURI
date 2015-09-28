#!/usr/bin/env python
from falcon_kit.FastaReader import FastaReader
import re
import sys
import os
import argparse

def parseHMMout(infilename, monomer_file_name, orientation):
  ID="NO ID"
  outfile=open(monomer_file_name, 'w')
  with open(infilename) as f:
    for line in f:
      l = line.strip().split()    
      if l == []:
        continue
      if ">>" in line:      
        ID = l[1]
      else:
        if "!" in line:
          low = int(l[12]) - 1
          high = int(l[13]) - 1
          if high-low+1 >= mono_len_threshold:
            outfile.write(">%s/%d_%d/%s\n" % (ID,low+1,high+1,orientation))
            outfile.write("%s\n" % seq_db[ID][low:high+1])
            #print ID,"\t",low+1,"\t",high+1,"\t",orientation
  outfile.close()

class DefaultList(list):
    def __copy__(self):
        return []

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('fasta_file',
                    action='store',
                    help ='Name of the FASTA file containing the reads')

parser.add_argument('hmm_file_fwd',
                    action = 'store',
                    help = 'Name of the forward HMM file')

parser.add_argument('hmm_file_rev',
                    action = 'store',
                    help = 'Name of the reverse complement HMM file')


parser.add_argument('-l',
                    action = 'store',
                    type = int,
                    default = 150,
                    dest = 'min_monomer_length',
                    help = 'Minimum monomer length')

parser.add_argument('--version', action='version', version='%(prog)s 0.2')

results = parser.parse_args()
in_seq_file = results.fasta_file
hmm_model_fwd = results.hmm_file_fwd
hmm_model_rev = results.hmm_file_rev
mono_len_threshold = results.min_monomer_length
monomers_file=in_seq_file.replace(".fa","_inferred_monomers.fa")

# Call hmmsearch, build hmms based on consensus alignments
os.system("rm -f hmmoutF.tbl hmmoutF.out; hmmsearch --cpu 8 --tblout hmmoutF.tbl -o hmmoutF.out  --notextw %s %s" % (hmm_model_fwd, in_seq_file)) 
os.system("rm -f hmmoutR.tbl hmmoutR.out; hmmsearch --cpu 8 --tblout hmmoutR.tbl -o hmmoutR.out  --notextw %s %s" % (hmm_model_rev, in_seq_file)) 

seq_db = {}
for r in FastaReader(in_seq_file):
    seq = r.sequence
    seq_db[r.name] = seq

parseHMMout("hmmoutF.out", "inferred_monomers_F.zzz", "F")
parseHMMout("hmmoutR.out", "inferred_monomers_R.zzz", "R")
os.system("cat inferred_monomers_F.zzz inferred_monomers_R.zzz > inferred_monomers.fa; rm inferred_monomers_F.zzz inferred_monomers_R.zzz")
