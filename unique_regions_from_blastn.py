#!/usr/bin/python

# This script is intended to identify unique regions of a specified
# minimum size from an all versus all blastn report. It wasn't written
# with general applicability in mind and has a number of idiosyncracies
# that make it non-trivial to generalize. One of these is that it relies
# on the lengths of the query sequences being encoded in the identifiers
# of the queries. This is the case for contigs produced by certain
# sequencing read assemblers such as the SPAdes assembler. The purpose
# of this script is to identify 1000 nt regions (hard-coded below)
# within a conting which is unique to that contig among all the contigs
# of the assembly. Starting from an assembly where the sequence
# identifiers satisify the abovementioned requirement for encoding
# length the next step is to run an all versus all blastn search. The
# script expects a non-canonical tab-delimited blastn output. For
# example:
# First make a BLAST database of your contigs
# $ makeblastdb -in contigs.fa -dbtype nucl -title contigs -out contigs
# Then run an all versus all blastn search
# $ blastn -task blastn -db contigs -query contigs.fa -evalue 1e-10 \
#   -out AVAblastn.txt -outfmt "6 qseqid sseqid pident length mismatch \
#   gapopen qstart qend sstart send evalue bitscore qlen slen gaps"
# 
# Then this script can be run on that output (hardcoded as AVAblastn.txt
# below)

import sys

class blastHit:
    def __init__(self, line, SEP='\t'):
        """Parses a custom output from NCBI's blast+ tool kit

        The output parsed is produced by running a blast search with the
        following output parameter:
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps'
        """
        line = line.rstrip('\n')
        cols = line.split(SEP)
        self.qseqid = cols[0]
        self.sseqid = cols[1]
        self.pident = float(cols[2])
        self.length = int(cols[3])
        self.mismatch = int(cols[4])
        self.gapopen = int(cols[5])
        self.qstart = int(cols[6])
        self.qend = int(cols[7])
        self.sstart = int(cols[8])
        self.send = int(cols[9])
        self.evalue = float(cols[10])
        self.bitscore = float(cols[11])
        self.qlen = int(cols[12])
        self.slen = int(cols[13])
        self.gaps = int(cols[14])

    def isSelfHit(self):
        return self.qseqid == self.sseqid

    def isFullSelfHit(self):
        if self.isSelfHit() and \
                self.qstart == 1 and \
                self.qend == self.qlen and \
                self.sstart == 1 and \
                self.send == self.slen and \
                self.qlen == self.slen and \
                self.mismatch == 0:
            return True
        else:
            return False

    def isFullReverse(self):
        if self.qlen == self.slen and \
                self.qstart == 1 and \
                self.qend == self.qlen and \
                self.sstart == self.slen and \
                self.send == 1 and \
                self.mismatch == 0:
            return True
        else:
            return False

    def areSameLength(self):
        """Tests if the subject and query are the same length

        This will probably usually need to be combined with testing that it
        isn't also a self hit.
        """
        return self.qlen == self.slen

if __name__ == '__main__':
    import re

    def sizeOfContig(contigName):
        return int(re.search(r'_length_(\d+)', contigName).group(1))
    
    def findContiguous(trueFalseList, minLen):
        """Returns 1-based coordinates of contiguous regions
        
        The trueFalseList is a list instance where each element is either True
        or False. The indices of the list correspond to positions on a
        continuous linear sequence. When a series of consecutive True values
        occur in a length greater than the minimum length, then the whole
        series will be returned in one-based coordinates. Separate contiguous
        regions will be comma-delimited
        """
        regions = []
        trueFalseList = [False] + trueFalseList
        i = 0
        while i < len(trueFalseList):
            if trueFalseList[i] == False:
                i += 1
            else:
                j = i
                while trueFalseList[j] == True:
                    j += 1
                    if j == len(trueFalseList):
                        break
                if j - i >= minLen:
                    regions.append([str(i), str(j-1)])
                i = j + 1
        return(regions)

    infile = 'AVAblastn.txt'
    lines = list(open(infile,'r').readlines())
    contigNames = list(set([line.split('\t')[0] for line in lines]))
    nrPositions = {contigName : [True for i in range(sizeOfContig(contigName))]
            for contigName in contigNames}
    for line in lines:
        bh = blastHit(line)
        if bh.isSelfHit():
            continue
        else:
            i1 = bh.qstart - 1
            i2 = bh.qend
            nrpos = nrPositions[bh.qseqid]
            for i in range(i1, i2):
                nrpos[i] = False
            nrPositions[bh.qseqid] = nrpos
    for contig in contigNames:
        nrpos = findContiguous(nrPositions[contig], 1000)
        if len(nrpos) == 0:
            print("%s None" % contig)
        else:
            print("%s" % (' '.join([contig + ':' + '-'.join(nr) for nr in nrpos])))


