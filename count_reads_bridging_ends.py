#! /usr/bin/python3

# Author: Patrick Denis Browne
# e-mail: pdbr@plen.ku.dk
# Licensed under the terms of the GNU General Public License v3.0

import re
maxInsertSize = 500

class sam_entry:
    def __init__(self, line):
        if line.startswith('@'):
            raise ValueError('Only headless SAM files are supported.')
        line = line.rstrip('\n')
        self.__cols = line.split('\t')
        assert len(self.__cols) >= 11, "Line has too few fields for SAM format\n%s" % (line)
        for i in [1, 3, 4, 7, 8]:
            self.__cols[i] = int(self.__cols[i])
        self.qname = self.__cols[0]
        self.flag = self.__cols[1]
        self.rname = self.__cols[2]
        self.pos = self.__cols[3]
        self.mapq = self.__cols[4]
        self.cigar = self.__cols[5]
        self.rnext = self.__cols[6]
        self.pnext = self.__cols[7]
        self.tlen = self.__cols[8]
        self.seq = self.__cols[9]
        self.qual = self.__cols[10]
    
    def __str__(self):
        return('\t'.join([str(col) for col in self.__cols]))

class cigar: # Cannot parse no CIGAR string represeted by '*'
    def __init__(self, cigarString):
        '''
        This parses the CIGAR strings present in SAM format alignments.
        If no CIGAR string was calculated for the alignment it will be
        represented by a '*'. This class will not parse such a string
        correctly.
        '''
        self.cigarString = cigarString
        if not self.cigarString == '*':
            self.subCigarStrings = re.findall(r'\d+[MINDSHPX=]', self.cigarString)
            self.parts = [[a[-1:], int(a[:-1])] for a in self.subCigarStrings]
            self.readLength = self.__readLength()
            self.alnLength = self.__alnLength()
        else:
            self.subCigarStrings = None
            self.parts = None
            self.readLength = None
            self.alnLength = None
    
    def __readLength(self):
        '''
        Calculates the length of the read based on the CIGAR string
        '''
        l = 0
        for part in self.parts:
            if part[0] in ['M', 'I', 'S', 'X', '=']:
                l += part[1]
        return(l)
    
    def __alnLength(self):
        '''
        Calculates the length of the alignment based on the CIGAR string
        '''
        l = 0
        for part in self.parts:
            if part[0] in ['M', 'D', 'N', '=', 'X']:
                l += part[1]
        return(l)
    
    def readMatchStartStop(self):
        '''
        Returns the minimum and maximum coordinates along the read at
        which are aligned to the reference.
        '''
        if self.readLength == None:
            return(None)
        poss = [i+1 for i in range(self.readLength)]
        i = 0
        subs = []
        for part in self.parts:
            if part[0] == 'S':
                for j in range(part[1]):
                    subs.append(i)
                    i += 1
            if part[0] in ['M', 'I', 'X', '=']:
                for j in range(part[1]):
                    i += 1
        subs.reverse()
        for i in subs:
            _ = poss.pop(i)
        if len(poss) == 0:
            return(None)
        return(min(poss), max(poss))
        
        

class sam_per_read:
    """
    Input:   A headless SAM file of namesorted reads. This could be made
             with samtools' sort function with the '-n' flag enabled. If
             any other filtering, e.g. based on bitwise flags or if you
             want all entries of each read mapped to a single template,
             the data needs to be preprocessed with samtools or another
             suitable tool. Python will generally be much too slow for
             this kind of task.
    
    Returns: An iterable yielding 2D lists of SAM-format alignments by
             read-name. The number of lists in the returned list will be
             the number of alignments for that read.
    """
    def __init__(self,filePath):
        self.fileh = open(filePath,'r')
        self.cursor_pos = self.fileh.tell()
    
    def __iter__(self):
        return(self)
    
    def next(self):
        currlines = [self.fileh.readline().rstrip('\n').split('\t')]
        self.cursor_pos = self.fileh.tell()
        if currlines == [['']]:
            raise StopIteration
        while self.fileh.readline().rstrip('\n').split('\t')[0] == currlines[0][0]:
            self.fileh.seek(self.cursor_pos)
            currlines.append(self.fileh.readline().rstrip('\n').split('\t'))
            self.cursor_pos = self.fileh.tell()
        self.fileh.seek(self.cursor_pos)
        return(currlines)
        
        

if __name__ == "__main__":
    import sys
    usage = '''reads_overlap_point.py aln.sam
    It only looks for reads overlapping the origin.
    The SAM file must be in plain text (BAM not supported) and headless
    (lines starting with '@' throw an exception).
    '''
    if '-h' in sys.argv or '--h' in sys.argv or '-help' in sys.argv or \
    '--help' in sys.argv:
        print(usage)
        quit()
    if len(sys.argv) > 2:
        print(usage)
        quit()
    
    def lenFromContigName(contigNameStr):
        a = re.search(r'_length_(\d+)_', contigNameStr)
        if a == None:
            return(None)
        else:
            return(int(a.group(1)))
    
    def isR1(sam_entry_obj):
        '''
        If the "segment unmapped" (0x4) and the "SEQ being reverse
        complemented" flags are NOT set, this returns True. Otherwise
        this returns False.
        '''
        if sam_entry_obj.flag & 0x4 == 0 and \
        sam_entry_obj.flag & 0x10 == 0:
            return(True)
        else:
            return(False)
    
    def isR2(sam_entry_obj):
        '''
        If the "segment unmapped" (0x4) flas is NOT set and the "SEQ
        being reverse complemented" flags are IS set, this returns True.
        Otherwise this returns False.
        '''
        if sam_entry_obj.flag & 0x4 == 0 and \
        sam_entry_obj.flag & 0x10 == 1:
            return(True)
        else:
            return(False)
    
    def sameReadForward(a, b):
        '''
        a and b are both sam_entry objects.
        It will be assessed if:
            a and b are both mapped
            a and b are both R1 or both R2
            a and b both map in the forward orientation
        '''
        if a.flag & 0x4 == 0 and b.flag & 0x4 == 0 and \
        a.flag & 0x40 == b.flag & 0x40 and \
        a.flag & 0x10 == 0 and b.flag & 0x10 == 0:
            return(True)
        else:
            return(False)
    
    def sameReadReverse(a, b):
        '''
        a and b are both sam_entry objects.
        It will be assessed if:
            a and b are both mapped
            a and b are both R1 or both R2
            a and b both map in the reverse orientation
        '''
        if a.flag & 0x4 == 0 and b.flag & 0x4 == 0 and \
        a.flag & 0x40 == b.flag & 0x40 and \
        a.flag & 0x10 != 0 and b.flag & 0x10 == 1:
            return(True)
        else:
            return(False)
    
    def pairedReads(a, b):
        '''
        a and b are both sam_entry objects.
        It will be assessed if:
            a and b are both mapped
            a and b are both paired
            (a is first segment and b is second segment) OR
            (a is second segment and b is first segment)
        It will return False if any of the above evaluate as False
        It will return 'R1R2' if a is R1 and b is R2.
        It will return 'R2R1' if a is R2 and b is R1.
        '''
        if a.flag & 0x1 != 0 and b.flag & 0x1 != 0 and \
        a.flag & 0x4 != 0 and b.flat & 0x4 != 0:
            pass
        else:
            return(False)
        if a.flag & 0x40 == 0 and b.flag & 0x80 == 0:
            return('R1R2')
        elif a.flag & 0x80 == 0 and b.flag & 0x40 == 0:
            return('R2R1')
        else:
            return(False)
    
    def singleReadSegmentBridges(a,b):
        '''
        a and b are both sam_entry objects
        '''
        # See if a and b are validly mapped, either both R1s or R2s and
        # mapped in the same orientat
        if sameReadForward(a, b) == True:
            strand = 'forward'
        elif sameReadReverse(a, b) == True:
            strand = 'reverse'
        else:
            return(False)
        
        # Make sure a and b are mapped to the same template
        if a.rname != b.rname:
            return(False)
        template_length = lenFromContigName(a.rname)
        
        aCigar = cigar(a.cigar)
        bCigar = cigar(b.cigar)
        read_length = aCigar.readLength
        
        if aCigar.readMatchStartStop() == None \
        or bCigar.readMatchStartStop() == None:
            return(False)
        
        if aCigar.readMatchStartStop()[1] < bCigar.readMatchStartStop()[0] and \
        a.pos > b.pos + bCigar.alnLength and strand == 'forward':
            if template_length - a.pos + b.pos + bCigar.alnLength <= read_length:
                return(True)
        
        if bCigar.readMatchStartStop()[1] < aCigar.readMatchStartStop()[0] and \
        b.pos > a.pos + aCigar.alnLength and strand == 'forward':
            if template_length - b.pos + a.pos + aCigar.alnLength <= read_length:
                return(True)
        
        if aCigar.readMatchStartStop()[1] < bCigar.readMatchStartStop()[0] and \
        a.pos + aCigar.alnLength < b.pos and strand == 'reverse':
            if template_length - b.pos + a.pos + aCigar.alnLength <= read_length:
                return(True)
        
        if bCigar.readMatchStartStop()[1] < aCigar.readMatchStartStop()[0] and \
        b.pos + bCigar.alnLength < a.pos and strand == 'reverse':
            if template_length - a.pos + b.pos + bCigar.alnLength <= read_length:
                return(True)
        
        return(False)
    
    def pairedReadBridge(a, b, maxInsertSize):
        if a.rname != b.rname:
            return(False)
        template_length = lenFromContigName(a.rname)
        
        aCigar = cigar(a.cigar)
        bCigar = cigar(b.cigar)
        aAlnLen = aCigar.readMatchStartStop()[1] - aCigar.readMatchStartStop()[0] + 1
        bAlnLen = bCigar.readMatchStartStop()[1] - bCigar.readMatchStartStop()[0] + 1
        
        if pairedReads(a, b) == 'R1R2':
            if a.flag & 0x10 != 0 and b.flag & 0x10 == 0 and \
            a.pos > b.pos + bAlnLen and \
            template_length - a.pos + b.pos + bAlnLen <= maxInsertSize:
                return(True)
            if a.flag & 0x10 == 0 & b.flag & 0x10 != 0 and \
            b.pos > a.pos + aAlnLen and \
            template_length - b.pos + a.pos + aAlnLen <= maxInsertSize:
                return(True)
        
        elif pairedReads(a, b) == 'R2R1':
            if b.flag & 0x10 != 0 and a.flag & 0x10 == 0 and \
            b.pos > a.pos + aAlnLen and \
            template_length - b.pos + a.pos + aAlnLen <= maxInsertSize:
                return(True)
            if b.flag & 0x10 == 0 & a.flag & 0x10 != 0 and \
            a.pos > b.pos + bAlnLen and \
            template_length - a.pos + b.pos + bAlnLen <= maxInsertSize:
                return(True)
        
        else:
            return(False)
        
    
    samf = sys.argv[1]
    samh = sam_per_read(samf)
    total = 0
    singleReadHits = 0
    pairedReadHits = 0
    supportedBySingleRead = {}
    supportedByPairedRead = {}
    while samh:
        try:
            readaln = samh.next()
            total += 1
            alns = [sam_entry('\t'.join(line)) for line in readaln]
            for i in range(len(alns)-1):
                for j in range(i, len(alns)):
                    if singleReadSegmentBridges(alns[i], alns[j]) == True:
                        print(alns[i])
                        print(alns[j])
                        print("\n")
                        supportedBySingleRead[alns[i].rname] = \
                        supportedBySingleRead.get(alns[i].rname, 0) + 1
                        singleReadHits += 1
                    if pairedReadBridge(alns[i], alns[j], maxInsertSize) == True:
                        print(alns[i])
                        print(alns[j])
                        print("\n")
                        supportedByPairedRead[alns[i].rname] = \
                        supportedByPairedRead.get(alns[i], 0) + 1
                        pairedReadHits += 1
        except StopIteration:
            break
    print("Supported by single reads")
    for k,v in supportedBySingleRead.items():
        print(k, v)
    print("\nSupported by paired reads")
    for k,v in supportedByPairedRead.items():
        print(k, v)
    print("\n")
    print("Total: %i\nSingleHits: %i\nPairedHits: %i" % (total, singleReadHits, pairedReadHits))
                    
    
