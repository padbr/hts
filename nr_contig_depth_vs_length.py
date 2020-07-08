#! /usr/bin/python

# This assess how deeply sequencing reads from a direct mobilome library
# and an MDA amplified mobilome library map to a mobilome assembly.
#
# It requires:
# (i) an assembled mobilome (this was specifically written for
# application to a metaplasmid SPAdes assembly of paired Illumina data).
# The relevant file is assigned to the 'CONTIGS' variable below.
# (ii) A file detailing the unique regions of the contigs of the
# assembly. Such a file could be obtained as detailed in
# https://github.com/padbr/hts/blob/master/unique_regions_from_blastn.py
# The format of the file should be evident in the example lines below
# where contigs with no unique regions are followed by 'None':
# Contig_1297_length_1635 None
# Contig_0318_length_2999:1-1751 Contig_0318_length_2999:1927-2999
# Contig_2544_length_1572 None
# Contig_0145_length_4769:1-4769
# The relevant file is assigned to the 'NR_REGIONS' variable below
# (iii) Mapping of two sequencing datasets to the assembly. The mapping
# could probably be done with any suitable program that maps fastq reads
# to a fasta assembly and produces SAM format output. Following mapping,
# the output SAM-format files should be processed with samtools to
# remove unmapped reads and supplementary alignments which can be
# achieved with `samtools view -F 0x804' and convert to BAM format.
# The relevant files are assigned to the 'DIRBAM' and 'AMPBAM' variables
# below. 

import re
import subprocess
import numpy as np

CONTIGS = 'contigs.fasta'
NR_REGIONS = 'nrRegions_min_1kb.txt'
DIRBAM = 'dir.sorted.aln.bam'
AMPBAM = 'amp.sorted.aln.bam'

def countReads(line, bamF):
    line = line.strip()
    a = line.split(' ')
    b = tuple(['samtools', 'view', bamF] + a)
    ps1 = subprocess.Popen(b, stdout=subprocess.PIPE)
    ps2 = subprocess.Popen(('cut', '-f', '1'), stdin=ps1.stdout, stdout=subprocess.PIPE)
    ps3 = subprocess.Popen(('sort'), stdin=ps2.stdout, stdout=subprocess.PIPE)
    ps4 = subprocess.Popen(('uniq'), stdin=ps3.stdout, stdout=subprocess.PIPE)
    try:
        strOutput = subprocess.check_output(('grep', '-c', '', '-'), stdin=ps4.stdout)
        ps1.wait()
        ps2.wait()
        ps3.wait()
        ps4.wait()
        output = int(strOutput.rstrip('\n'))
        return output
    except:
        return int(0)

def bamDepth(line, bamF, trim=0):
    """
    line: A Region recognized by samtools. It must specify start and stop
          positions, regardless of whether or not the whole region should
          be matched
    bamF: A bam file. It must be sorted and indexed
    trim: To avoid edge effects in mapping, it is probably a good idea to
          trim 50 to 100 nucleotides off the edge of the regions to only
          look at positions with more robust coverage values.

    Returns:
    length:The number of positions for which coverage is being reported
    average: The average coverage depth
    median: The median coverage depth
    """
    line = line.rstrip('\n')
    regions = line.split(' ')
    depths = []
    numPositions = 0
    for region in regions:
        cName, dStart, dStop = re.search(r'^([^:]+):(\d+)-(\d+)$', region).groups()
        assert int(dStart) <= int(dStop), "The start wasn't before the stop\n%s" % (region)
        dStart = str(int(dStart) + trim)
        dStop = str(int(dStop) - trim)
        assert int(dStart) <= int(dStop), "The trimming was too much\n%s" % (region)
        newRegion = cName + ':' + dStart + '-' + dStop
        matchLen = 1 + int(dStop) - int(dStart)
        numPositions += matchLen
        ps1 = subprocess.Popen(('samtools', 'view', '-h', bamF, newRegion), stdout=subprocess.PIPE)
        output = subprocess.check_output(('samtools', 'depth', '-'), stdin=ps1.stdout)
        ps1.wait()
        output = output.rstrip('\n')
        outlines = output.split('\n')
        if outlines == ['']:
            outdepths = [0]
        else:
            outdepths = [int(outline.split('\t')[2]) for outline in outlines]
        depths = depths + outdepths
    assert numPositions > 0, "Depth cannot be calculated for this region (zero or negative length)\n%s" % (line)
    rMean = float(sum(depths)) / numPositions
    rMedian = np.median(depths + [0]*(numPositions - len(depths)))
    return numPositions, rMean, rMedian

def lenFromRecName(recname):
    rlen = int(re.search(r'.*_length_(\d+)', recname).group(1))
    return rlen

# rlens = {str(rec.id) : len(rec) for rec in SeqIO.parse(CONTIGS, 'fasta')}

print('\t'.join(['Name', 'Length', 'Dir Num Reads', 'Amp Num Reads', 'Dir Avg depth', 'Amp Avg depth', 'Dir Median depth', 'Amp Median depth', 'NR length']))
for line in open(NR_REGIONS, 'r').readlines():
    line = line.rstrip('\n')
    if line.endswith('None'):
        continue
    recName = line.split(':')[0]
    recLen = lenFromRecName(recName)
    dirnumReads = countReads(line, DIRBAM)
    matchLen, dirmeanDepth, dirmedianDepth = bamDepth(line, DIRBAM, trim=50)
    ampnumReads = countReads(line, AMPBAM)
    matchLen, ampmeanDepth, ampmedianDepth = bamDepth(line, AMPBAM, trim=50)
    print("%s\t%i\t%i\t%i\t%s\t%s\t%s\t%s\t%i" % (recName, recLen, dirnumReads, ampnumReads, str(dirmeanDepth), str(ampmeanDepth), str(dirmedianDepth), str(ampmedianDepth), matchLen))


