#!/usr/bin/python

'''Report the longest common substring between two fasta files. We will calculate the longest string for the forward and reverse of the first file compared to the forward of the second file'''

import sys
sys.path.append('/home3/redwards/bioinformatics/Modules')
import rob

try:
    file1 = sys.argv[1]
    file2 = sys.argv[2]
except:
    sys.stderr.write(sys.argv[0] + " <file 1> <file 2>\n")
    sys.exit(-1)


def longestCommonSubstring(s1, s2):
    '''This is taken straight from the wikibooks page, and is creating a matrix to look up. Dynamic programming'''
    m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in xrange(1, 1 + len(s1)):
        for y in xrange(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]

fa1=rob.readFasta(file1)
fa2=rob.readFasta(file2)

longest = ""
for id1 in fa1.keys():
    for id2 in fa2.keys(): 
        sys.stderr.write("Comparing " + id1 + " to " + id2 + "\n")
        test = longestCommonSubstring(fa1[id1], fa2[id2])
        if len(test) > len(longest):
            longest = test
        sys.stderr.write("Comparing rc " + id1 + " to " + id2 + "\n")
        test = longestCommonSubstring(rob.rc(fa1[id1]), fa2[id2])
        if len(test) > len(longest):
            longest = test

print "\t".join([file1, file2, len(longest), longest])

