from falcon_kit import kup, falcon, DWA, get_consensus, get_alignment
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sys
import pickle
from pbcore.io import FastaReader


print "RUNNING IN DEBUG MODE: Adjust j range!"
arguments = len(sys.argv) - 1
if arguments != 2:
    print "Need two arguments: Start and end"    
else:
    start = int(sys.argv[1])
    end = int(sys.argv[2])

if start >= end:
    print "start > end\n"
    quit()

monomer_list_hupac = []
monomer_list_chm = []

for r in FastaReader("pread_HuPac_monomers.fa"):
    monomer_list_hupac.append(r.sequence)    
uniq_hupac = set(monomer_list_hupac)
       
for r in FastaReader("pread_CHM1_10x_monomers.fa"):
   monomer_list_chm.append(r.sequence)    
uniq_chm = set(monomer_list_chm)


onlyinH = uniq_hupac.difference(uniq_chm)
onlyinC = uniq_chm.difference(uniq_hupac)
common = uniq_chm.intersection(uniq_hupac)

labelH = ["H" for i in range(len(onlyinH))]
labelC = ["C" for i in range(len(onlyinC))]
labelBoth = ["B" for i in range(len(common))]

labeledOnlyH = zip(labelH, onlyinH)
labeledOnlyC = zip(labelC, onlyinC)
labelBoth = zip(labelBoth, common)

alllabeled = labeledOnlyH + labeledOnlyC + labelBoth

aln_data = []

#!!!
labeli, monoi = alllabeled[1]
labelj, monoj = alllabeled[2]
alignment = DWA.align(monoi, len(monoi), monoj, len(monoj), 50, 1)      
size = alignment[0].aln_str_size
dist = alignment[0].dist        
aln_str1 = alignment[0].q_aln_str
aln_str0 = alignment[0].t_aln_str                
#!!!

for i in range(start, end): #range(len(alllabeled)):
    print i
    for j in range(i+1, 200): #range(i+1, len(alllabeled)): 
        labeli, monoi = alllabeled[i]
        labelj, monoj = alllabeled[j]
        #!!!alignment = DWA.align(monoi, len(monoi), monoj, len(monoj), 50, 1)        
        #!!!size = alignment[0].aln_str_size
        #!!!dist = alignment[0].dist        
        #!!!aln_str1 = alignment[0].q_aln_str
        #!!!aln_str0 = alignment[0].t_aln_str                
        
        #print aln_str0, "\n", aln_str1
        #print "size = ", size, "\n dist = ", dist
        
        if size != 0:
            score = 1.0*dist/size
            aln_data.append( (i, j, round(score,3)) )
            
        #!!!DWA.free_alignment(alignment)            

filename = "alignment_scores_" + str(start).zfill(6) + "-" + str(end).zfill(6) + ".pickle"
pfile = open(filename, "wb+")
pickle.dump(aln_data, pfile)
pfile.close()

for item in aln_data:
    print item
