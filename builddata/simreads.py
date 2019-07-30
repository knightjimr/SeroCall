
import sys
import re
if len(sys.argv) != 2:
    sys.stderr.write("Error:  symreads.py fastaFile\n")
    sys.exit(-1)

lastacc = ""
lastseq = ""

#mode = "Newbler"
mode = "BWA"

if mode == "Newbler":
    outfp = open("simreads_altorder.fastq", "w")
else:
    outfp1 = open("simreads_R1.fastq", "w")
    outfp2 = open("simreads_R2.fastq", "w")

step = 1

fp = open(sys.argv[1], "r")
while True:
    line = fp.readline()

    if line.startswith(">") or not line:
        seqlen = len(lastseq)
        if seqlen > 0:
            for pos in range(1, seqlen - 200, step):
                fseq = lastseq[pos:pos+100]

                if mode == "Newbler":
                    outfp.write("@%s_%d_%d/1\n%s\n+\n%s\n" % (lastacc, pos, pos+200, fseq, "A" * 100))
                else:
                    outfp1.write("@%s_%d_%d\n%s\n+\n%s\n" % (lastacc, pos, pos+200, fseq, "A" * 100))

                rseq = ""
                for x in reversed(range(pos+101,pos+201)):
                    ch = lastseq[x]
                    if ch == 'A':
                        ch = 'T'
                    elif ch == 'C':
                        ch = 'G'
                    elif ch == 'G':
                        ch = 'C'
                    elif ch == 'T':
                        ch = 'A'
                    elif ch == 'N':
                        ch = 'N'
                    else:
                        sys.stderr.write("Error:  Invalid seq char:  %s\n", ch)
                        sys.exit(-1)
                    rseq += ch

                if mode == "Newbler":
                    outfp.write("@%s_%d_%d/2\n%s\n+\n%s\n" % (lastacc, pos, pos+200, rseq, "A" * 100))
                else:
                    outfp2.write("@%s_%d_%d\n%s\n+\n%s\n" % (lastacc, pos, pos+200, rseq, "A" * 100))

        if not line:
            break

        match = re.match(">([^\s]+)\s", line)
        lastacc = match.group(1)
        lastseq = ""
        sys.stderr.write("Reading %s...\n" % lastacc)
        sys.stderr.flush()
    else:
        lastseq += line.strip()
fp.close()
