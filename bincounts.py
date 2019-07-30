
import sys
import re

# bincounts.py
#
# This is an internal script to parse the lines of the bam file, and compute the read counts occurring
# in each bin of the references (and difference locations).  It is not intended for use as a command-line program.
#

def main():
    name = "sample"
    capsuleFile = ""
    otherFile = ""
    arg = 1
    while arg < len(sys.argv):
        if arg + 1 < len(sys.argv) and sys.argv[arg] == "-n":
            name = sys.argv[arg+1]
            arg += 2
        elif arg + 1 < len(sys.argv) and sys.argv[arg] == "-c":
            capsuleFile = sys.argv[arg+1]
            arg += 2
        elif arg + 1 < len(sys.argv) and sys.argv[arg] == "-o":
            otherFile = sys.argv[arg+1]
            arg += 2
        else:
            break

    if len(sys.argv) - arg != 2:
        sys.stderr.write("Usage:  bincounts [-n sampleName] [-c capsule.sam] [-o other.sam] serobins.txt { align.sam | - }\n")
        sys.exit(-1)

    bins = []
    totalcounts = []
    uniqcounts = []

    binmap = {}
    pileups = {}
    pileupBins = {}

    # Read the "serobins" bin design, and setup the data structures for holding reads counts for each bin.
    # Also, setup pileups for any pseudogene difference location.

    fp = open(sys.argv[arg])
    for line in fp:
        f = line.rstrip("\n").split("\t")

        # Setup bins, totalcounts and uniqcounts
        
        bins.append((f[1], int(f[2]), int(f[3]), f[4] if len(f) > 4 else "", f[0]))
        totalcounts.append(0)
        uniqcounts.append(0)

        if f[1] not in binmap:
            binmap[f[1]] = []
        binmap[f[1]].append(len(bins) - 1)

        # Create any pileups, if a pseudogene.

        if f[0].startswith("pg"):
            if f[0] not in pileups:
                pileups[f[0]] = Pileup(bins[-1][0], bins[-1][1], bins[-1][2])
                pileupBins[f[0]] = []
            pileupBins[f[0]].append(len(bins)-1)
    fp.close()

    # Setup input read counts, and open the SAM file.

    numReads = 0
    numUnmapped = 0
    numGenome = 0
    numCapsule = 0
    numOther = 0
    
    if sys.argv[arg+1] == "-":
        fp = sys.stdin
    else:
        fp = open(sys.argv[arg+1])

    # Debug writing of unusual read alignments.

    capfp = open(capsuleFile, "w") if capsuleFile else None
    otherfp = open(otherFile, "w") if otherFile else None

    lastacc = ""
    lastf = []

    errcnts = [ 0 ] * 10

    # Read and process the SAM lines.

    for (cnt, line) in enumerate(fp):
        if line.startswith("@"):
            continue

        if cnt > 0 and cnt % 1000000 == 0:
            sys.stderr.write("  -> %d\n" % cnt)
            sys.stderr.flush()

        f = line.strip().split("\t")

        # Skip supplemental alignments.
        flags = int(f[1])
        if (flags & 0x800) or (flags & 0x100):
            continue

        # Collect the two primary alignments for the same readpair, generating an error if they don't exist.

        if not lastf:
            lastf = f
            continue

        if f[0] != lastf[0]:
            sys.stderr.write("Error:  Cound not find two primary alignments for read:  %s\n" % lastf[0])
            sys.exit(-1)

        f1 = lastf
        f2 = f

        lastf = None
        numReads += 1

        # Count unmapped reads, genomic matches or other unusual alignments.
         
        if (int(f1[1]) & 0x4) and (int(f2[1]) & 0x4):
            numUnmapped += 1
            continue

        if f1[2].startswith("Spneumo") and f2[2].startswith("Spneumo"):
            numGenome += 1
            continue

        if f1[2] == "*" or f1[5] == "*" or f1[6] != "=" or f1[7] == f1[3]:
            numOther += 1
            #if otherfp:
            #    otherfp.write("\t".join(f1) + "\n" + "\t".join(f2) + "\n")
            continue

        # For well-aligned readpairs, get the error counts for the reads, 
        # and check if above the threshold.

        errs = getErrCnt(f1) + getErrCnt(f2)
        if errs >= 10:
            numOther += 1
            if otherfp:
                otherfp.write("\t".join(f1) + "\n" + "\t".join(f2) + "\n")
            continue

        numCapsule += 1
        if capfp:
            capfp.write("\t".join(f1) + "\n" + "\t".join(f2) + "\n")

        # Build the histogram of per-readpair error counts.

        errcnts[errs] += 1

        # For each of the primary reads, add the read to the bin counts.

        for x in [ f1, f2 ]:
            serotype = x[2]
            pos = int(x[3])
            mapq = int(x[4])

            # Use the cigar string to compute the end pos for the alignment.

            cigar = x[5]
            endpos = pos
            while cigar:
                match = re.match("(\d+)([HDIMS])", cigar)
                if match is None:
                    sys.stderr.write("Error:  Invalid cigar string:  %s\n" % x[5])
                    sys.exit(-1)

                if match.group(2) in ("M", "D"):
                    endpos += int(match.group(1))
                cigar = cigar[match.end():]

            # For each bin for the aligned serotype, add the read count to each overlapping bin.

            if serotype not in binmap:
                sys.stderr.write("Error:  Serotype not in list:  %s\n" % serotype)
                sys.exit(-1)

            for i in binmap[serotype]:
                if not (pos < bins[i][2] and bins[i][1] < endpos):
                    continue
                
                # For the small difference locations, ensure that the overlap is a complete overlap,
                # to avoid edge effects from the alignment counts.

                if bins[i][2] - bins[i][1] < 50 and not (pos <= bins[i][1] and bins[i][2] <= endpos):
                    continue

                # If the bin is not a pseudogene difference location, just increment the readcounts.

                if not bins[i][4].startswith("pg"):
                    totalcounts[i] += 1
                    if mapq > 0:
                        uniqcounts[i] += 1
                    continue
                    
                # If the bin is a pseudogene difference location, build the alignment and add to the pileup

                align = BamAlign(x)
                if bins[i][3]:
                    adjustAlign(align, bins[i][3])
                pileups[bins[i][4]].add(align)

    fp.close()

    # For each pseudogene pileup, call the pseudogene percentages.

    for pg in pileups:
        (fcnt, pcnt) = pileups[pg].callPseudogene()
        for bidx in pileupBins[pg]:
            if bins[bidx][3]:
                totalcounts[bidx] = pcnt
                uniqcounts[bidx] = pcnt
            else:
                totalcounts[bidx] = fcnt
                uniqcounts[bidx] = fcnt


    # Generate the output.

    l = [ "NumReads", "NumUnmapped", "NumGenome", "NumOther", "NumCapsule" ]
    l2 = [ numReads, numUnmapped, numGenome, numOther, numCapsule ]

    sys.stdout.write("##fileformat=BinCountsv1.0\n")
    for i in range(len(l)):
        sys.stdout.write("##%s=%d\n" % (l[i], l2[i]))
    sys.stdout.write("#SAMPLE\tSEROTYPE\tSTART\tEND\tTOTAL\tUNIQUE\n")
    for i in range(len(bins)):
        sys.stdout.write("%s\t%s\t%d\t%d\t%d\t%d%s\n" % (name, bins[i][0], bins[i][1], bins[i][2], totalcounts[i], uniqcounts[i],
                                                         "\t%s" % bins[i][4] if bins[i][4].startswith("diff") or bins[i][4].startswith("pg") else ""))
                            

def getErrCnt(fields):
    cnt = 0

    # Count any hard or soft clipped bases as errors.

    cigar = fields[5]
    match = re.match("(\d+)[HS]", cigar)
    if match is not None:
        cnt += int(match.group(1))

    match = re.search("(\d+)[HS]$", cigar)
    if match is not None:
        cnt += int(match.group(1))

    # Add the error count from the alignment itself.

    for x in fields:
        if x.startswith("NM:i:"):
            cnt += int(x[5:])

    return cnt

# adjustAlign
#
# For pseudogene alignments where the read aligned to the serotype with the pseudogene (i.e., the serotype has 
# the variant that causes loss-of-function of the gene), readjust the alignment to be the alignment of the
# functional version of the gene (i.e., without the variant in the reference, but now with the variant in the
# query).  This way, even though the reads may align to any of the serotype sequences, the pileup works from
# the functional gene's sequence.
def adjustAlign(align, diffstr):
    match = re.match(":(\d+)(.+)>(.+)", diffstr)
    diffpos = int(match.group(1))
    qstr = match.group(2)
    rstr = match.group(3)

    rpos = align.refstartpos
    i = 0
    while i < len(align.qalign):
        if rpos == diffpos:
            if qstr.startswith("-"):
                if align.qalign[i] != '-':
                    align.qalign = align.qalign[:i] + [ch for ch in qstr ] + align.qalign[i:]
                    align.ralign = align.ralign[:i] + [ch for ch in rstr ] + align.ralign[i:]
                    align.refendpos += 1
            elif rstr.startswith("-"):
                if align.ralign[i] != '-':
                    align.qalign = align.qalign[:i] + [ch for ch in qstr ] + align.qalign[i:]
                    align.ralign = align.ralign[:i] + [ch for ch in rstr ] + align.ralign[i:]
                    align.refendpos += 1
            else:
                align.ralign[i] = rstr

            break

        if align.ralign != '-':
            rpos += 1
        i += 1

# Class to handle the storing of pseudogene pileup statistics.

class Pileup:
    def __init__(self, accno, start, end):
        self.accno = accno
        self.start = start
        self.end = end

        reglen = end - start
        self.refch = [ 'X' ] * reglen
        self.depth = [ 0 ] * reglen
        self.snpcnt = []
        self.indelcnt = []
        for i in range(reglen):
            self.snpcnt.append([ 0 ] * 5)
            self.indelcnt.append({})

    # Add a new alignment to the pileup counts.

    def add(self, align):
        nuc2idx = { 'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3, 'N' : 4 }

        rpos = align.refstartpos
        i = 0
        while i < len(align.qalign):
            if align.qalign[i] == '-':
                j = i + 1
                while j < len(align.qalign) and align.qalign[j] == '-':
                    j += 1

                altstr = '-' + "".join(align.ralign[i:j])
                if rpos >= self.start and rpos < self.end:
                    if altstr not in self.indelcnt[rpos-self.start]:
                        self.indelcnt[rpos-self.start][altstr] = 0
                    self.indelcnt[rpos-self.start][altstr] += 1

                while i < j:
                    if rpos >= self.start and rpos < self.end:
                        self.depth[rpos-self.start] += 1
                    rpos += 1
                    i += 1
            elif align.ralign[i] == '-':
                j = i + 1
                while j < len(align.ralign) and align.ralign[j] == '-':
                    j += 1

                altstr = '+' + "".join(align.qalign[i:j])
                if rpos >= self.start and rpos < self.end:
                    if altstr not in self.indelcnt[rpos-self.start]:
                        self.indelcnt[rpos-self.start][altstr] = 0
                    self.indelcnt[rpos-self.start][altstr] += 1

                i = j
            else:
                if rpos >= self.start and rpos < self.end:
                    if self.refch[rpos-self.start] == 'X':
                        self.refch[rpos-self.start] = align.ralign[i]
                    self.depth[rpos-self.start] += 1
                    self.snpcnt[rpos-self.start][nuc2idx[align.qalign[i]]] += 1

                rpos += 1
                i += 1

    # Compute the proportions of functional gene vs. non-functional gene from the pileup counts,
    # looking for any frameshift insertion or stop gained SNP in the non-reference counts.

    def callPseudogene(self):
        idx2nuc = [ 'A', 'C', 'G', 'T', 'N' ]

        maxlof = 0.0
        maxdepth = 0
        for i in range(len(self.refch)):

            # Check any SNP that is not the reference, and has a readcount >= 4.

            for j in range(len(self.snpcnt[i])):
                if self.snpcnt[i][j] >= 4 and idx2nuc[j] != self.refch[i]:
                    #sys.stderr.write("%d:  %s>%s  %.1f (%d of %d)\n" % (self.start+i, self.refch[i], idx2nuc[j], self.snpcnt[i][j] * 1.0 / self.depth[i], self.snpcnt[i][j], self.depth[i]))

                    # Check the variant codon to see if it is a stop codon.

                    codonidx = i - (i % 3)
                    if codonidx + 2 < len(self.refch):
                        codon = (self.refch[codonidx] if codonidx != i else idx2nuc[j]) + \
                                (self.refch[codonidx+1] if codonidx+1 != i else idx2nuc[j]) + \
                                (self.refch[codonidx+2] if codonidx+2 != i else idx2nuc[j])
                        if codon == "TGA":
                            lof = self.snpcnt[i][j] * 1.0 / self.depth[i]
                            if lof > maxlof:
                                maxlof = lof
                                maxdepth = self.depth[i]

            # Check any frameshifting variant that has a readcount >= 4.

            for indel in self.indelcnt[i]:
                if self.indelcnt[i][indel] >= 4:
                    #sys.stderr.write("%d:  %s>%s  %.1f (%d of %d)\n" % (self.start+i, self.refch[i], indel, self.indelcnt[i][indel] * 100.0 / self.depth[i], self.indelcnt[i][indel], self.depth[i]))
                    if (len(indel)-1) % 3 != 0:
                        lof = self.indelcnt[i][indel] * 1.0 / self.depth[i]
                        if lof > maxlof:
                            maxlof = lof
                            maxdepth = self.depth[i]

        if maxdepth == 0:
            maxdepth = int(sum(self.depth) / len(self.depth))

        # Compute the read counts of functional and pseudogene using the max depth of any LOF variant.

        fcnt = int(maxdepth * (1.0 - maxlof))
        pcnt = int(maxdepth * maxlof)

        return (fcnt, pcnt)

# Class to handle parsing of the SAM alignments, producing alignment strings for the query and reference.

class BamAlign:
    def __init__(self, f):
        self.alignFlag = False

        flags = int(f[1])

        self.suppFlag = ((flags & 0x900) != 0)
        self.strand = "<" if flags & 0x10 else ">"

        self.alignFlag = True
        self.refchr = f[2]
        self.mq = int(f[4])
        
        rpos = int(f[3])
        qseq = f[9]
        qscoreseq = f[10]

        mdtag = ""
        for i in range(11,len(f)):
            if f[i].startswith("MD:Z:"):
                mdtag = f[i][5:]

        ralign = []
        qalign = []
        qscore = []

        prefix = ""
        suffix = ""
        prefFlag = True

        self.refstartpos = rpos

        spos = 0 
        cigar = f[5]
        while cigar:
            match = re.match("(\d+)([MDHIS])", cigar)
            if match is None:
                sys.stderr.write("Error:  Invalid cigar string:  %s\n%s\n" % (f[5], "\t".join(f)))
                sys.exit(-1)

            num = int(match.group(1))
            ch = match.group(2)

            if ch == 'M':
                prefFlag = False
                for x in range(num):
                    qalign.append(qseq[spos])
                    qscore.append(qscoreseq[spos])
                    ralign.append(qseq[spos])
                    rpos += 1
                    spos += 1
            elif ch == 'D':
                prefFlag = False
                for x in range(num):
                    qalign.append("-")
                    qscore.append("!")
                    ralign.append("N")
                    rpos += 1
            elif ch == 'I':
                prefFlag = False
                for x in range(num):
                    qalign.append(qseq[spos])
                    qscore.append(qscoreseq[spos])
                    spos += 1
                    ralign.append("-")
            elif ch == "S":
                if prefFlag:
                    prefix += "%dS" % num
                else:
                    suffix += "%dS" % num
                spos += num
            elif ch == "H":
                if prefFlag:
                    prefix += "%dH" % num
                else:
                    suffix += "%dH" % num
                pass
            else:
                sys.stderr.write("Error:  Invalid cigar string:  %s\n%s\n" % (f[5], "\t".join(f)))
                sys.exit(-1)

            cigar = cigar[match.end():]

        self.refendpos = rpos - 1

        if mdtag:
            i = 0
            apos = 0
            while i < len(mdtag) and apos < len(ralign):
                num = 0
                while i < len(mdtag) and mdtag[i].isdigit():
                    num = num * 10 + int(mdtag[i])
                    i += 1

                if i == len(mdtag):
                    break

                x = 0
                while x < num and apos < len(ralign):
                    if ralign[apos] == "-":
                        apos += 1
                    else:
                        x += 1
                        apos += 1
                while apos < len(ralign) and ralign[apos] == "-":
                    apos += 1
                if apos == len(ralign):
                    break

                if mdtag[i] != '^':
                    ralign[apos] = mdtag[i]
                    apos += 1
                    i += 1
                else:
                    i += 1
                    while i < len(mdtag) and apos < len(ralign) and not mdtag[i].isdigit():
                        if qalign[apos] != "-":
                            sys.stderr.write("Error:  Unable to parse MD string, pos %d, apos %d:  %s\n%s\n   %s\n   %s\n" % (i, apos, mdtag, "\t".join(f), qalign, ralign))
                            sys.exit(-1)

                        ralign[apos] = mdtag[i]
                        apos += 1
                        i += 1

        self.qalign = qalign
        self.ralign = ralign
        self.qscore = qscore
        self.prefix = prefix
        self.suffix = suffix

if __name__ == "__main__":
    main()
