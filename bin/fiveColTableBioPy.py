import sys,os,re
from Bio import SeqIO

try:
    gb_file=sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <genbank file>")


gb_record = SeqIO.read(open(gb_file,"r"), "genbank")
print ">Feature " + gb_record.name
for f in gb_record.features:
    if f.type == "source":
        # I am not printing the source as it has a bunch of errors
        continue
    beg = f.location.start.real + 1
    end = f.location.end.real
    if f.location.strand == -1:
        beg, end = end, beg
    print "\t".join([str(beg), str(end), f.type])
    for q in f.qualifiers:
        if f.type == "tRNA" and q == "product":
            pr = re.sub('-[ATGC]{3}', '', f.qualifiers[q][0])
            pr=pr.rstrip()
            print "\t\t\t" + str(q) + "\t" + pr
        elif q == "product":
            # remove the EC number from the product name
            pr = re.sub('\(EC\s+[\d\.\-]+\)', '', f.qualifiers[q][0])
            pr=pr.rstrip()
            print "\t\t\t" + str(q) + "\t" + pr
        else:
            print "\t\t\t" + str(q) + "\t" + str(f.qualifiers[q][0])

