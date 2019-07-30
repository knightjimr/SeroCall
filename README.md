# SeroCall
SeroCall can identify and quantitate the capsular serotypes in Illumina whole-genome sequencing samples of *S. pneumoniae*,
calculating abundances of each serotype in mixed cultures.  The software is written in Python (compatible with Python 2 or 3), and is freely available under an open source GPLv3 license.

Note (July 30, 2019):  The software has just been released, and is still considered in beta testing.

## Contents
  * [Installation and Dependencies](#installation)
  * [Running SeroCall](#running)
  * [Building/Updating the Reference Database](#updating)
  * [Citation](#citation)
------------------------------

## Installation and Dependencies

SeroCall has the following dependencies:
* Python 2.7+ or 3+ (tested using 2.7.11 and 3.5.1)
* BWA version >= 0.7.15
* Samtools version >= 1.3 (only used if you change/update the reference database files)

Currently, to install SeroCall, download or clone this github repository, then add the SeroCall directory to your
PATH environment variable.  Installation through GitHub releases, pip, bioconda and Docker are forthcoming.


## Running SeroCall

```
usage: serocall [-h] [-o OUTPUT] [-t THREAD] r1file r2file

positional arguments:
  r1file                R1 FASTQ file (can be gzipped)
  r2file                R2 FASTQ file (can be gzipped)

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Prefix for output files [default 'sero']
  -t THREAD, --thread THREAD
                        Number of processors to use [default 1]
```

In its simplest form, the serocall command just takes the the R1 and R2 fastq files from a paired-end sequencing run.
The files can be gzipped or uncompressed, and there are no limitations or patterns to the given filenames.

A "-t #" option can be used to have BWA MEM use multiple threads for processing (BWA MEM is very efficient with
parallelizing the alignments, so using multiple threads will reduce the running time of the tool).

SeroCall will generate two output files, a *_counts.txt containing the intermediate "bin" counts, counting how many reads
aligned across the regions of the serotype references, and a *_calls.txt containing the final calls.  By default, the
output file names will use the prefix "sero" (outputting files "sero_counts.txt" and "sero_calls.txt"), but that can
be changed using the "-o prefix" option.

The following is an example sero_calls.txt file:
```
##fileformat=SeroCallv1.0
##NumReads=3420303
##NumUnmapped=129037
##NumGenome=3232094
##NumOther=31298
##NumCapsule=27874
#SEROTYPE       PERCENTAGE
19F     61.8%
01      38.2%
```
The first line of the output is a file format tag.  The other lines, before the "#SEROTYPE..." header, contain metrics about
the read data, including the following:
* NumReads - the number of input reads
* NumUnmapped - number of reads not mapped to serotype or *S. pneumoniae* genomes
* NumGenome - number of reads mapped to the *S. pneumoniae* genomes (not the serotype sequences)
* NumOther - number of other alignments, including partial, high error and chimeric alignments
* NumCapsule - number of reads well-aligned to the serotype capsular sequences (i.e., used to call serotypes)

The "#SEROTYPE..." header and following lines are tab-delimited lines reporting the serotypes called and their percentages.
If no lines appear after the header, no serotypes were called from the data.

The sero_counts.txt file is technically an intermediate file, but can be useful for post-analysis QC of the calls.  The
initial lines of an example file are the following:
```
##fileformat=BinCountsv1.0
##NumReads=3420303
##NumUnmapped=129037
##NumGenome=3232094
##NumOther=31298
##NumCapsule=27874
#SAMPLE SEROTYPE        START   END     TOTAL   UNIQUE
sero    01      1       500     373     369
sero    01      501     1000    546     508
sero    01      1001    1500    612     584
sero    01      1501    2000    640     339
sero    01      2001    2500    512     322
sero    01      2501    3000    556     556
sero    01      3001    3500    568     551
sero    01      3501    4000    669     669
sero    01      4001    4500    637     637
sero    01      4501    5000    852     852
sero    01      5001    5500    775     775
sero    01      5501    6000    698     698
sero    01      6001    6500    655     655
sero    01      6501    7000    729     729
sero    01      7001    7500    786     786
sero    01      7501    8000    732     732
sero    01      8001    8500    763     763
sero    01      8501    9000    780     780
sero    01      9001    9500    729     729
.
.
.
```

After the header lines, the rest of the file is a tab-delimited file identifying the regions of the serotype sequences,
and then reporting the total counts of aligned reads, and the count of uniquely aligning reads in each region.  The final
lines in the file contain an extra column with a "diff" identifier, marking them as difference locations for
distinguishing closely related serogroups.

## Building/Updating the Reference Database

*WARNING:  This functionality has not been tested externally, does require some (to a lot of) care when making changes, and
should be needed by only a small subset of users.*

The current software has pre-built databases containing the same serotype sequences as the PneumoCaT v1.2 database, and
building the database is not required in order to use the serocall command.  The package does provide the ability
to change and update the references used by SeroCall, but this has not been extensively tested.

The builddata sub-directory contains the data files and scripts necessary to rebuild the SeroCall database files (found
in the "data" sub-directory of the repository).  In that sub-directory, the serotype capsular sequences are found in the
"serotypes" sub-directory, the S. pneumoniae references are found in the "genomes" sub-directory and the difference
locations are found in the "pcat_diffs.txt" file.

To add a new serotype capsular sequence to the database, create a FASTA file containing only that sequence, and then place
it in the serotypes sub-directory.  IMPORTANT:  The file must end with ".fa" and the name of the file must match the
accession given on the FASTA header.  So, for a new serotype like "10X", the name of the file should be "10X.fa" and the
first line of the file should be ">10X ..." where after the first space on the line, any text is allowed.

To add a new genome reference to the database, create a masked FASTA file of the reference sequence and place it in
the genomes sub-directory.  Any name can be given to the file, however the accession on the FASTA header must begin with
"Spneumo_".  IMPORTANT:  The capsular sequence must be masked from the genome sequence (or all of the reads will align
to it, confusing SeroCall).  To do that, megablast/blat/minimap2 the serotype sequences against the genome reference,
identify any region of the genome with very high identity (>98%) longer than 100 bp, and replace those regions' bases
with N's in the FASTA file.  (Remember, the goal is not to match to this sequence, but to use it as a decoy so that
genomic reads don't align to the serotype sequences.)

If you believe you want to make changes/additions to the difference locations, please email j.knight@yale.edu for
instructions on how to do that (there are a number of subtle issues that must be accounted for, in order for the software
to properly identify serotypes using these difference).

Once the changes have been made to these files, then the command "bash builddata.sh" command will rebuild the database
files and install them into the top-level "data" directory, for use by the serocall command.  Things to note when
running this command:
* Both bwa and samtools are required for this step, so they should be in your PATH.
* All of the initial construction of the files occurs in the "build" directory, then a final transfer occurs into the
top-level "data" directory.
* A backup of the existing top-level "data" directory files is made, placed in "data/bak".  (But, only a single backup
is made, so if builddata.sh is run multiple times, the previous backup is removed.)
* This command will download a pair of FASTQ files, totalling 1.5 GB, and store them in the "genomeOnly" sub-directory
the first time it is run.  These files contain a real run of a capsule-free sample, and is used to identify those
serotype regions where genomic reads may align into the serotypes.

## Citation

To cite the SeroCall package, please use the following citation:

   Knight J, Dunne E, Mulholland E, Saha S, Satzke C, Tothpal A, Weinberger D.  Determining the Serotype Composition of Mixed Samples of Pneumococcus using Whole Genome Sequencing.  To be submitted.

