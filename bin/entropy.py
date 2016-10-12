import matplotlib
matplotlib.use("AGG")
from rob import *
import matplotlib.pyplot as mpl
import sys

if len(sys.argv) < 2:
    sys.exit(sys.argv[0] + " <kmer size> <fasta file>")

shan = {}

k = int(sys.argv[1])
for file in sys.argv[2:]:
    sys.stderr.write(file + "\n")
    fa = readFasta(file)
    shan[file]={}
    count=0
    for idd in fa:
        count = count+1
        #if count > 10000:
        #    break
        s = "%0.2f" % (shannon(fa[idd], k))
        if s in shan[file]:
            shan[file][s] = shan[file][s]+1
        else:
            shan[file][s] = 1


colors = ['g', 'r', 'c', 'm', 'y', 'k']
for file in shan:
    c=colors.pop()
    x=[]
    y=[]
    for d in shan[file]:
        x.append(d)
        y.append(shan[file][d])
    labelname = file.split('/')
    label=labelname[len(labelname)-1]
    mpl.plot(x,y, c, ls='', marker='.', label=label)

mpl.title("Entropy with " + str(k) + "-mers")
mpl.legend(loc=2)
mpl.savefig('entropy.' + str(k) + '.png')
