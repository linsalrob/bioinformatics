"""

Translate a sequence (eg. a metagenome) and print the longest ORF
"""


import sys
import os
from robseq import translate
from rob import read_fasta
from rob import rc

try:
    faf = sys.argv[1]
    size = int(sys.argv[2])
except:
    sys.exit(sys.argv[0] + " <fasta file> <min orf size in amino acids>")



fa = read_fasta(faf)

for seqid in fa:
    outputid = seqid.split(' ')[0]
    c = 0
    orfs = []
    location = {}
    for frame in range(3):
        dna = fa[seqid][frame:]
        prot = translate(dna)
        pieces = prot.split('*')
        orfs += pieces
        for p in pieces:
            location[p] = " frame: " + str(frame) + " strand: +"

    original = rc(fa[seqid])
    for frame in range(3):
        dna = original[frame:]
        prot = translate(dna)
        pieces = prot.split('*')
        orfs += pieces
        for p in pieces:
            location[p] = " frame: -" + str(frame) + " strand: -"

    longest = max(orfs, key=len)
    if len(longest) > size:
        print(">" + outputid + " " + location[longest] + "\n" + longest)

