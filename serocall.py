
import sys
import argparse
from itertools import groupby
from math import log

# Process the command-line arguments

parser = argparse.ArgumentParser(prog="serocall")
parser.add_argument("-o", "--output", type=str, default="sero", help="Output file prefix")
parser.add_argument("-t", "--threshold", type=float, default=1.0, help="Calling threshold %")
parser.add_argument("-d", "--debug", type=str, default="", help="Debug serotype list")
parser.add_argument("expected", type=str, help="File of expected counts")
parser.add_argument("observed", type=str, help="File of observed counts")

args = parser.parse_args()

### Read the baseline file, get the serotypes and bins, and then construct the baseline and genome baseline counts.

# Read the baseline file, dividing the lines into genome and serotype lines.

bLines = []
gLines = []
fp = open(args.expected)
for line in fp:
    if line.startswith("#"): continue

    f = line.rstrip("\n").split("\t")
    if f[0] == "genome":
        gLines.append(f)
    else:
        bLines.append(f)
fp.close()

# From the genome line, extract the serotypes and the bins. Each is an ordered list, but with map lookups
 
serotypes = []
seromap = {}

bins = []
binmap = {}

for f in gLines:
    serotype = f[1]
    startpos = f[2]
    endpos = f[3]

    if serotype not in seromap:
        seromap[serotype] = len(serotypes)
        serotypes.append(serotype)

    sidx = seromap[serotype]

    key = ":::".join(f[1:4])   # Bin lookup is the tuple (serotype, startpos, endpos)
    binmap[key] = len(bins)
    bins.append((sidx, startpos, endpos))

numbins = len(bins)

debug_serolist = []
if args.debug:
    l = args.debug.split(",")
    for s in l:
        sname = s.upper()
        if sname not in seromap:
            sys.stderr.write("Error:  Invalid serotype for debugging:  %s\n" % s)
            sys.exit(-1)
        debug_serolist.append(seromap[sname])

# Extract the genome baseline from the genome lines.

genomeBaseline = [ (int(f[4]), int(f[5])) for f in gLines ]

# Construct the serotype baselines from the serotype lines.

baseline = []
maxER = 0
maxEU = 0

for f in bLines:
    serotype = f[0]
    try:
        binER = int(f[4])
        binEU = int(f[5])
    except:
        sys.stderr.write("%s\n" % str(f))
        raise

    key = ":::".join(f[1:4])
    bidx = binmap[key]
    sidx = seromap[serotype]

    if binER > maxER: maxER = binER
    if binEU > maxEU: maxEU = binEU

    baseline.append((sidx, bidx, binER, binEU, f[6] if len(f) == 7 else ""))

#### Read in the observed counts, as well as the metrics lines and sample name.

sampleName = ""
metricsLines = []

OR = []
OU = []

fp = open(args.observed)
bidx = 0
for line in fp:
    if line.startswith("#"):
        if not (line.startswith("#SAMPLE") or line.startswith("##fileformat=")):
            metricsLines.append(line)
        continue

    f = line.rstrip("\n").split("\t")

    if not sampleName:
        sampleName = f[0]

    key = ":::".join(f[1:4])
    if key not in binmap or binmap[key] != bidx:
        sys.stderr.write("Error:  Baseline and sample bin counts differ.\n")
        sys.exit(-1)

    OR.append(int(f[4]))
    OU.append(int(f[5]))
    bidx += 1
fp.close()

#### Initialize the factor levels, then run the simulated annealing rounds.

F = [ 1.0 ] * len(serotypes)
MR = [ 0.0 ] * numbins
MU = [ 0.0 ] * numbins
RR = [ 0.0 ] * numbins
RU = [ 0.0 ] * numbins

for round in range(1, 101):
    #if round % 10 == 0:
    #    sys.stdout.write("Round %d...\n" % round)

    # Given the serotype factor levels, compute the expected mixture counts per bin using the expected counts.
    for i in range(numbins):
        MR[i] = 0.0
        MU[i] = 0.0

    for (sidx, bidx, ER, EU, bintype) in baseline:
        if F[sidx] < 0.0001: continue

        MR[bidx] += ER * F[sidx]
        MU[bidx] += EU * F[sidx]

    # Compute the observed / expected ratios for each bin.
    for i in range(numbins):
        RR[i] = OR[i] / MR[i] if MR[i] > 0 else 0.0
        RU[i] = OU[i] / MU[i] if MU[i] > 0 else 0.0

    # Compute weighted totals/cnts of ratios, based on per-bin percentage unique reads and total reads.
    # Also compute the ratio of unique mapping reads to total reads.

    TR = []
    TU = []
    for i in range(len(serotypes)):
        TR.append([])
        TU.append([])

    Uratio = [ 0.0 ] * len(serotypes)
    Uratio_cnt = [ 0 ] * len(serotypes)

    for (sidx, bidx, ER, EU, snpGroup) in baseline:
        # Skip cross-serotype baseline entries, genome positive locations and SNP regions
        if sidx != bins[bidx][0]: continue
        if genomeBaseline[bidx][0] >= 10: continue
        if snpGroup: continue

        if F[sidx] < 0.0001: continue
    
        #if sidx in debug_serolist:
        #    print bidx, RR[bidx], RU[bidx], OR[bidx], OU[bidx], ER, EU, genomeBaseline[bidx]

        # Compute weighted mean of ratios, with weights as percent ER/EU counts from maximum
        wr = ER * 1.0 / maxER
        wu = EU * 1.0 / maxEU

        if ER > 40:
            TR[sidx].append((min(10.0, RR[bidx]), wr))
        if EU > 40:
            TU[sidx].append((min(10.0, RU[bidx]), wu))

        # For the R/U ratio, compute the ratio of unique reads to total mapped reads, then weight that by the
        # ratio of total reads (for the bin) over the max total reads (i.e., how many there should have been if
        # this bin was all uniquely mapping reads).

        wratio = EU * 1.0 / ER
        Uratio[sidx] += wr * wratio
        Uratio_cnt[sidx] += wr

    # Adjust the abundances based on the new ratios.

    Fprev = list(F)
    for sidx in range(len(serotypes)):
        if F[sidx] < 0.0001:
            continue

        tr = TR[sidx]
        tu = TU[sidx]

        if len(tr) < 4:
            ratioR = -1.0
        else:
            tr.sort()

            if len(tr) >= 8:
                winsize = int(len(tr) * 0.6) if len(tr) >= 8 else len(tr)
            elif len(tr) >= 5:
                winsize = len(tr) - 2
            elif len(tr) >= 3:
                winsize = len(tr) - 1
            else:
                winsize = len(tr)

            mini = 0
            mindiff = tr[winsize-1][0] - tr[0][0]
            for i in range(1, len(tr)):
                if i + winsize - 1 >= len(tr):
                    break

                d = tr[i+winsize-1][0] - tr[i][0]
                if d < mindiff:
                    mini = i
                    mindiff = d

            minir = mini

            cnt = sum([ tr[i][1] for i in range(mini, mini+winsize) ])
            ratioR = -1.0
            c = 0.0
            for i in range(mini, mini+winsize):
                c += tr[i][1]
                if c * 2.0 >= cnt:
                    ratioR = tr[i][0]
                    break

        if len(tu) < 4:
            ratioU = -1.0
        else:
            tu.sort()

            if len(tu) >= 8:
                winsize = int(len(tu) * 0.6) if len(tu) >= 8 else len(tu)
            elif len(tu) >= 5:
                winsize = len(tu) - 2
            elif len(tu) >= 3:
                winsize = len(tu) - 1
            else:
                winsize = len(tu)

            wutotal = sum([ t[1] for t in tu ])
            wzero = sum([ t[1] for t in tu if t[0] == 0.0 ])

            start = 0 if wzero * 2 > wutotal else next(i for i in range(len(tu)) if tu[i][0] != 0.0)

            if start + winsize > len(tu):
                winsize = len(tu) - start

            mini = start
            mindiff = tu[start+winsize-1][0] - tr[start][0]
            for i in range(start+1, len(tu)):
                if i + winsize - 1 >= len(tu):
                    break

                d = tu[i+winsize-1][0] - tu[i][0]
                if d < mindiff:
                    mini = i
                    mindiff = d

            miniu = mini

            cnt = sum([ tu[i][1] for i in range(mini, mini+winsize) ])
            ratioU = -1.0
            c = 0.0
            for i in range(mini, mini+winsize):
                c += tu[i][1]
                if c * 2.0 >= cnt:
                    ratioU = tu[i][0]
                    break

        if ratioR < 0.0 and ratioU < 0.0:
            UvsRratio = 0.5
            finalRatio = 1.0
        elif ratioR < 0.0:
            UvsRratio = 1.0
            finalRatio = ratioU
        elif ratioU < 0.0:
            UvsRatio = 0.0
            finalRatio = ratioR
        else:
            UvsRratio = max([ t[1] for t in tu ])  # Uratio[sidx] / Uratio_cnt[sidx] if Uratio_cnt[sidx] > 0 else 0.0
            finalRatio = ratioR * (1.0 - UvsRratio) + ratioU * UvsRratio

        origF = F[sidx]
        newF = F[sidx] * finalRatio
        F[sidx] += (newF - F[sidx]) * 0.5  # Only make half the adjustment in this round.
        if newF < 0.0002 or F[sidx] < 0.002:
            F[sidx] = 0.0

        if sidx in debug_serolist:
            #print tr
            print tu
            cnt = 0.0
            for (i, t) in enumerate(tr):
                print i, t, cnt
                cnt += t[1]
            print "Round %d, %s - ratioR=%.5f  ratioU=%.5f   UvsR=%.4f  finalRatio=%.5f  %.5f -> %.5f (%.5f)" % (round, serotypes[sidx], ratioR, ratioU, UvsRratio, finalRatio, origF, F[sidx], newF)
            print "--------------------------------------------------"

    # Stop if nothing has changed.

    #if debug_serolist:
    #    print F
    #    print "--------------------------------------------------"
        
    numChanged = sum([ 1 for i in range(len(serotypes)) if abs(Fprev[i] - F[i]) > 0.000001 ])
    
    if numChanged == 0:
        break

#### For the serogroups that require SNP groups to distinguish between them, readjust the factor levels
#### based on the SNP group counts.

snpGroupBins = [ (snpGroup, sidx, bidx, ER, EU) for (sidx, bidx, ER, EU, snpGroup) in baseline if snpGroup and sidx == bins[bidx][0] ]
snpGroupBins.sort()

snpGroupSerotypes = set()
Fsnp = [ 1000000.0 ] * len(serotypes)

for key, groupiter in groupby(snpGroupBins, lambda x: x[0]):
    group = list(groupiter)
    s = serotypes[group[0][1]]
    if s in ("06A", "06B", "06C", "06D", "11A", "11D", "11B", "11C"):
        l = list(group)
    else:
        l = [ t for t in group if F[t[1]] != 0.0 ] # list(group)
        if len(l) == 1: continue

    # Skip if none of the serotypes in the group have a positive factor level. 
    if sum([ F[t[1]] for t in l ]) == 0:
        continue

    # If this is a new serogroup, record the set of serotypes in the group as a string.
    serotypesInGroup = ":::".join([ str(t[1]) for t in l ])
    if serotypesInGroup not in snpGroupSerotypes:
        snpGroupSerotypes.add(serotypesInGroup)

    # Compute the ratio of observed over expected for each serotype in the SNP group,
    # then keep the minimal ratio for each serotype.
    for (snpGroup, sidx, bidx, ER, EU) in l:
        if ER == 0:
            continue

        print snpGroup, serotypes[sidx], bidx, ER, EU, OR[bidx], OU[bidx]

        ratio = 0.0
        if EU > 0:
            ratio = OU[bidx] * 1.0 / EU
        elif ER > 0:
            ratio = OR[bidx] * 1.0 / ER

        if ratio < Fsnp[sidx]:
            Fsnp[sidx] = ratio

ambiguousCalls = []
for serotypesInGroup in snpGroupSerotypes:
    serolist = [ int(x) for x in serotypesInGroup.split(":::") ]

    # Readjust the factor levels of the serotypes in the serogroup by the proportion of their
    # Fsnp ratios (as computed above).

    total = sum([ Fsnp[sidx] for sidx in serolist ])

    # If none of the SNP differences provided enough information to distinguish between
    # serotypes, then create an ambigous serotype (containing all members of the group) as the call.
    if total < 0.0001:
        ambName = "/".join([ serotypes[sidx].replace("serotype_", "") for sidx in serolist ])
        Famb = sum([ F[sidx] for sidx in serolist ])
        ambiguousCalls.append([ serolist[0], Famb, ambName ])
        for sidx in serolist:
            F[sidx] = 0.0
        continue

    # Convert the Fsnp values for the serotypes into percentages, then readjust the factor levels as 
    # those percentages of the total factor level for the group.

    pcts = [ Fsnp[sidx] / total for sidx in serolist ]
    Ftotal = sum([ F[sidx] for sidx in serolist ])

    for i in range(len(serolist)):
        sidx = serolist[i]
        F[sidx] = pcts[i] * Ftotal

#### Make 15B/15C, 24A/24F and 32A/32F ambiguous calls, to match PneumoCaT reporting

pcatAmbiguousCalls = [ [ '15B', '15C' ], [ '24A', '24B', '24F' ], [ '32A', '32F' ] ]
for s in pcatAmbiguousCalls:
    serolist = [ seromap[serotype] for serotype in s ]
    ambName = "/".join([ serotypes[sidx].replace("serotype_", "") for sidx in serolist ])
    Famb = sum([ F[sidx] for sidx in serolist ])
    if Famb < 0.002:
        continue
    ambiguousCalls.append([ serolist[0], Famb, ambName ])
    for sidx in serolist:
        F[sidx] = 0.0

#### Make the final calls, based on the final factor levels, plus any ambiguous calls.

allcalls = [ [ sidx, F[sidx], serotypes[sidx] ] for sidx in range(len(serotypes)) if F[sidx] > 0.002 ]
allcalls += ambiguousCalls

if len(allcalls) == 0:
    finalcalls = []
else:
    alltotal = sum([ t[1] for t in allcalls ])

    finalcalls = [ [ -t[1] ] + t for t in allcalls if t[1] * 100.0 / alltotal > args.threshold ]
    finaltotal = sum([ t[2] for t in finalcalls ])

finalcalls.sort()

outfp = open("%s_calls.txt" % args.output, "w")
outfp.write("##fileformat=SeroCallv1.0\n")
for line in metricsLines:
    outfp.write("%s" % line)
outfp.write("#SEROTYPE\tPERCENTAGE\n")
#if len(finalcalls) == 0:
#    outfp.write("No serotypes found.\n")
for (tmp, sidx, Fcall, serotype) in finalcalls:
    outfp.write("%s\t%.1f%%\n" % (serotype, Fcall * 100.0 / finaltotal))
    sys.stdout.write("%s\t%.1f%%\n" % (serotype, Fcall * 100.0 / finaltotal))
outfp.close()


