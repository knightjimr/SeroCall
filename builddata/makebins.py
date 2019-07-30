
import sys
import re

if len(sys.argv) != 4:
    sys.stderr.write("Usage:  makebins.py binSize serotypes.fasta.fai pcat_diffs.txt\n")
    sys.exit(-1)

binSize = int(sys.argv[1])

fp = open(sys.argv[2])
for line in fp:
    f = line.rstrip("\n").split("\t")
    stype = f[0]
    seqlen = int(f[1])

    pos = 0
    while pos + binSize < seqlen:
        sys.stdout.write("seq\t%s\t%d\t%d\n" % (stype, pos+1, pos+binSize))
        pos += binSize
    if pos + (binSize / 2) < seqlen:
        sys.stdout.write("seq\t%s\t%d\t%d\n" % (stype, pos+1, seqlen))
fp.close()

diffcnt = 0
fp = open(sys.argv[3])
for line in fp:
    if line.startswith("#"): continue

    f = line.rstrip("\n").split("\t")
    diffcnt += 1

    num = int(f[0])
    f = f[1:]

    pgflag = False
    ranges = []
    for i in range(num, num + num):
        m1 = re.match("p(\d+)\-(\d+)(.*)", f[i])
        if m1 is not None:
            ranges.append((int(m1.group(1)), int(m1.group(2)), m1.group(3)))
            pgflag = True
            continue

        m1 = re.match("(\d+)\-(\d+)", f[i])
        if m1 is not None:
            ranges.append((int(m1.group(1)), int(m1.group(2))))
            continue

        m2 = re.match("(\d+)", f[i])
        if m2 is None:
            sys.stderr.write("Error:  Invalid position column %d:  %s\n" % (i, str(f)))
            sys.exit(-1)

        pos = int(m2.group(1))
        ranges.append((pos - 10, pos + 10))

    for j in range(num):
        sys.stdout.write("%s%d\t%s\t%d\t%d%s\n" % ("diff" if not pgflag else "pg", diffcnt, f[j], ranges[j][0], ranges[j][1],
                                                   "" if not pgflag else "\t" + ranges[j][2]))
fp.close()
