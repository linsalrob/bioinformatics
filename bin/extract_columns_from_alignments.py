
from Bio import AlignIO
import sys

if len(sys.argv) < 2:
    sys.exit("alignment file to parse?")

filename = sys.argv[1]
align = AlignIO.read(filename, 'clustal')
seqname = filename.replace(".aln", "")

l = align.get_alignment_length()
ids = [record.id for record in align]
different = {}
for id in ids:
    different[id]={}

allInteresting=[]

for i in range(l):
    col = align[:,i]
    # how many letters are there in this column
    letters={}
    for l in col:
        letters[l]=1
    # ignore columns that have a single letter
    if len(letters.keys()) == 1:
        continue

    interestingPosn = i+1 # convert to 1-indexed string
    for idx in range(len(col)):
        different[ids[idx]][interestingPosn] = col[idx]
    allInteresting.append(interestingPosn)

print "Strain\t" + "\t".join([seqname + "_" + str(i) for i in allInteresting])
for id in ids:
    sys.stdout.write(id)
    for inter in allInteresting:
        sys.stdout.write("\t" + different[id][inter])
    print






