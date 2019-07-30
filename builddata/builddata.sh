
set -e

if ! [ -e genomeOnly/T4J_genomeOnly_R1_001.fastq.gz ] ; then
    cd genomeOnly
    wget 'http://sysg1.cs.yale.edu:3010/tW6omLeDs96s3t2R7ghvL31bw2yib/share/T4J_genomeOnly_*.fastq.gz'
    cd ..
fi 

rm -rf build tmp
mkdir -p build
mkdir -p tmp

cat serotypes/*.fa > tmp/serotypes.fasta
samtools faidx tmp/serotypes.fasta

cat tmp/serotypes.fasta genomes/* > build/serotypes_plus_multi_genome.fasta
samtools faidx build/serotypes_plus_multi_genome.fasta
bwa index build/serotypes_plus_multi_genome.fasta

python makebins.py 500 tmp/serotypes.fasta.fai pcat_diffs.txt > build/serobins_500.txt

for file in serotypes/*.fa ; do
  NAME=`head -1 $file | sed 's/>//'`
  echo $NAME
  cd tmp
  python ../simreads.py ../$file
  bwa mem -t 10 -I 200 -M ../build/serotypes_plus_multi_genome.fasta simreads_R1.fastq simreads_R2.fastq > bwa.sam
  python ../../bincounts.py -n $NAME ../build/serobins_500.txt bwa.sam >> serotype_baseline.txt
  rm -f simreads_R1.fastq simreads_R2.fastq bwa.sam
  cd ..
done

echo genomeOnly
cd tmp

bwa mem -t 10 -M ../build/serotypes_plus_multi_genome.fasta ../genomeOnly/T4J_genomeOnly* > bwa.sam
python ../../bincounts.py -n genome ../build/serobins_500.txt bwa.sam > genomeOnly_baseline.txt
rm -f bwa.sam

cat genomeOnly_baseline.txt serotype_baseline.txt > ../build/baseline_bin500.txt

cd ..

rm -rf ../data/bak
mkdir -p ../data/bak
mv ../data/*.* ../data/bak
mv build/* ../data

