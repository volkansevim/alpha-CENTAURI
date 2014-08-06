from falcon_kit import kup, falcon, DWA, get_consensus, get_alignment
import matplotlib.pyplot as plt
from pbcore.io import FastaReader

monomer_list_hupac = []
monomer_list_chm = []

for r in FastaReader("pread_HuPac_monomers.fa"):
    monomer_list_hupac.append(r.sequence)    
uniq_hupac = set(monomer_list_hupac)
       
for r in FastaReader("pread_CHM1_10x_monomers.fa"):
   monomer_list_chm.append(r.sequence)    
uniq_chm = set(monomer_list_chm)

print "\t\t\t HuPac\t CHM1_10x"
print "--------------------------------------"
print "Unique Monomer Count\t", len(uniq_hupac), "\t", len(uniq_chm)
print "Total Monomer Count\t", len(monomer_list_hupac), "\t", len(monomer_list_chm)
print "% Unique Monomers\t", 1.*len(uniq_hupac)/len(monomer_list_hupac), "\t", 1.*len(uniq_chm)/len(monomer_list_chm)

nonoverlapping = uniq_chm.symmetric_difference(uniq_hupac)
intersect = uniq_chm.intersection(uniq_hupac)
union = uniq_chm.union(uniq_hupac)
print "\nIntersection = ", len(intersect)
print "Union = ", len(union)
print "Nonoverlapping = ", len(nonoverlapping)
print "%Nonoverlapping = ", 1.*len(nonoverlapping)/len(union)
print

#

#print "alignment"

#onlyinH = uniq_hupac.difference(uniq_chm)
#onlyinC = uniq_chm.difference(uniq_hupac)
#common = uniq_chm.intersection(uniq_hupac)

#labelH = ["H" for i in range(len(onlyinH))]
#labelC = ["C" for i in range(len(onlyinC))]
#labelBoth = ["B" for i in range(len(common))]

#labeledOnlyH = zip(labelH, onlyinH)
#labeledOnlyC = zip(labelC, onlyinC)
#labelBoth = zip(labelBoth, common)

#alllabeled = labeledOnlyC + labeledOnlyH + labelBoth

#aln_data = []
#for i in range(10): #range(len(alllabeled)):
    #for j in range(10): #range(len(alllabeled)):
        #labeli, monoi = alllabeled[i]
        #labelj, monoj = alllabeled[j]
        #alignment = DWA.align(monoi, len(monoi), monoj, len(monoj), 50, 1)        
        #size = alignment[0].aln_str_size
        #dist = alignment[0].dist        
        #aln_str1 = alignment[0].q_aln_str
        #aln_str0 = alignment[0].t_aln_str                
        
        #print aln_str0, "\n", aln_str1
        #print "size = ", size, "\n dist = " dist
        
        #if size != 0:
            #score = 1.0*dist/size
            #aln_data.append( (i, j, score) )
            
#DWA.free_alignment()            

#print "alllabeled:", len(alllabeled)        
#print "aln_data: ", len(aln_data)
