import sys
import os
from robseq import translate
from rob import read_fasta

try:
    faf = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " fasta file")



fa = read_fasta(faf)

for i in fa:
    p = translate(fa[i])
    # print(i + "\t" + p)
    print(">" + i + "\n" + p)
