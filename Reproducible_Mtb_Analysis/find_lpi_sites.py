#!/usr/bin/env python3
from Bio import SeqIO
import re

motif = "[CG]G.TA.C[GC]"
nstart = 3 #Number of bases to extend the circular record at the start
nend = 3 #Number of bases to extend the circular record at the "end"
noffset = 4 #Number of bases to offset the match by (if a T is the last char and A is first char in the record then 0 is the match_index)
out_file = 'LPI_sites.csv'
fp = open(out_file, 'w')
fp.write('contig,insertion_site\n')
for r in SeqIO.parse("H37Rv_GenBank.gb", "genbank"):
    #print(repr(r.seq))
    extend_r = str(r.seq[-nstart:] + r.seq + r.seq[1:nend]) #Add nucleotides at "start" and "end" for circular chromosomes
    contig = r.id
    #contig = ""

    matches = re.finditer(motif,extend_r)
    for m in matches:
        match_index = m.start() - nstart + noffset
        fp.write(contig + ',' + str(match_index) + '\n')
        
