

import rob
import sys

# 1404927386.fasta  analyzed_sequences.txt  annotations.txt
#

faf=None

try:
    faf=sys.argv[1]
except IndexError:
    sys.stderr.write("Please provide a fasta file\n")
    sys.exit(0)


fa = rob.readFasta(faf)


analyzed=[]
with open('analyzed_sequences.txt', 'r') as asf:
    for line in asf:
        pieces=line.rstrip()
        analyzed.append(pieces)
        if pieces not in fa:
            sys.stderr.write(pieces + " has been analyzed but is not in " + faf + "\n")

for f in fa:
    if f not in analyzed:
        sys.stderr.write("NOT ANALYZED: " + f + "\n")


annotated=[]
with open('annotations.txt', 'r') as asf:
    for line in asf:
        pieces=line.split("\t")
        annotated.append(pieces[0])
        if pieces[0] not in fa:
            sys.stderr.write(pieces[0] + " has been annotated but is not in " + faf + "\n")

for f in fa:
    if f not in annotated:
        sys.stderr.write("NOT ANNOTATED: " + f + "\n")


