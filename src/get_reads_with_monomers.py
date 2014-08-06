#!/usr/bin/env python
from pbcore.io import FastaReader
import subprocess, shlex, os

fr = FastaReader("/mnt/secondary/Share/vsevim/centromere/HuPacData/0-fasta_files/000a5794_q.fa")
for record in fr:
	print record.name, len(record.sequence), record.md5

fr.close()
