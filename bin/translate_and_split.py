
import sys
import os
from robseq import translate
from rob import read_fasta
from rob import rc

try:
    faf = sys.argv[1]
    size = int(sys.argv[2])
except:
    sys.exit(sys.argv[0] + " <fasta file> <min orf size>")



fa = read_fasta(faf)

for seqid in fa:
    outputid = seqid.split(' ')[0]
    c = 0
    for frame in range(3):
        dna = fa[seqid][frame:]
        prot = translate(dna)
        pieces = prot.split('*')
        for p in pieces:
            if len(p) >= size:
                c +=1
                print(">" + outputid + "_" + str(c) + " frame: " + str(frame) + " strand: +\n" + p)


    original = rc(fa[seqid])
    for frame in range(3):
        dna = original[frame:]
        prot = translate(dna)
        pieces = prot.split('*')
        for p in pieces:
            if len(p) >= size:
                c +=1
                print(">" + outputid + "_" + str(c) + " frame: " + str(frame) + " strand: -\n" + p)


