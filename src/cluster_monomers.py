from falcon_kit import kup, falcon, DWA, get_consensus, get_alignment
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from pbcore.io import FastaReader

threshold = 0.1
G = nx.Graph()

for i, j, idt in aln_data:
    if idt < threshold:
        G.add_edge(i, j)
print "Total # clusters ", nx.number_connected_components(G)
clusterIdx = 0
clusterContent = []

for C in nx.connected_components(G):
    print "Component ", clusterIdx
    dummy = [0, 0 ,0]
    for idx in C:
        if alllabeled[idx][0] == "H":
            dummy[0] += 1
        elif alllabeled[idx][0] == "C":
            dummy[1] += 1
        elif alllabeled[idx][0] == "B":
            dummy[2] += 1
    clusterContent.append(dummy)        
    clusterIdx += 1    

f = open("output.txt", "w+")       
f.write("Cluster#\tOnly HuPac\tOnly CHM1\tBoth")
for i in range(len(clusterContent)):    
    f.write(str(i) + "\t" + str(clusterContent[i][0])+ "\t" + str(clusterContent[i][1])+ "\t" + str(clusterContent[i][2]) + "\n")
f.close()            
        
