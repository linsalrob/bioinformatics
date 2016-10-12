import sys
import rob
import robseq
from alignment import score_alignment, gapped_alignment

try:
    fname = sys.argv[1]
except:
    sys.exit(sys.argv[0] + " <fasta file>")

prot = 'ALMEGQTFDKSAYPKLAVAYPSGVIPDMRGQTIKGKPSGRAVLSAEADGVKAHSHSASASSTDLGTKTTSSFDYGTKGTNSTGGHTHSGSGSTSTNGEHSHYIEAWNGTGVGGNKMSSYAISYRAGGSNTNAAGNHSHTFSFGTSSAGDHSHSVGIGAHTHTVAIGSHGHTITVNSTGNTENTVKNIAFNYIVRLA'

fa = rob.read_fasta(fname)

#print(prot + "\n")
for s in fa:
    maxd = -99999
    best={}
    t = fa[s]
    for i in range(0,3):
        p = robseq.translate(t[i:])
        d = score_alignment(p, prot)
        #print(">" + s + " [edit distance: " + str(d) + "] [frame: " + str(i) + "]\n" + str(p))
        if d  > maxd:
            maxd = d
            best = {'frame' : i, 'seq' : p}
        elif d == maxd:
            sys.stderr.write("Two equal scores for " + s + "\n")
    t = rob.rc(fa[s])
    for i in range(0,3):
        p = robseq.translate(t[i:])
        d = score_alignment(p, prot)
        #print(">" + s + " [edit distance: " + str(d) + "] [frame: " + str(-i) + "]\n" + str(p))

        if d > maxd:
            maxd = d
            best = {'frame' : -i, 'seq' : p}
        elif d == maxd:
            sys.stderr.write("Two equal scores for " + s + "\n")
    #print("\n\n")

    (d, p1, p2) = gapped_alignment(prot, best['seq'])
    print(">" + s + " [edit distance: " + str(maxd) + "] [frame: " + str(best['frame']) + "]\n" + str(p1) + "\n" + str(p2) + "\n")

