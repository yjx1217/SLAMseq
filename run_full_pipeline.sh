#!/bin/bash
######### Set path for dependencies ##########
ngm_dir="/public/shared/tools/NextGenMap-0.5.5/bin/ngm-0.5.5"
samblaster_dir="/public/jxyue/Projects/SLAMseq/build/samblaster-v.0.1.26"
trimmomatic_dir="/public/jxyue/Projects/SLAMseq/build/Trimmomatic-0.38"
samtools_dir="/public/jxyue/Projects/SLAMseq/build/samtools-1.9"
htslib_dir="/public/jxyue/Projects/SLAMseq/build/samtools-1.9/htslib-1.9/"
picard_dir="/public/jxyue/Projects/SLAMseq/build/Picard-v2.19.0"
slamdunk_dir="/public/jxyue/Projects/SLAMseq/build/miniconda3/bin/"
varscan_dir="/public/jxyue/Projects/SLAMseq/build/varscan-2.4.5"
seqtk_dir="/public/jxyue/Projects/SLAMseq/build/seqtk-1.3"
java_dir="/usr/bin"
#############################################
PATH=$samtools_dir:$PATH

##### Set input file information ############
refseq="ref.genome.fa" # reference genome in FASTA format
sample_id="infected_4h" # sample id prefix for output files
R1_reads="SRR13573515_1.fastq.gz" # The R1 read file in FASTQ format 
R2_reads="SRR13573515_2.fastq.gz" # The R2 read file in FASTQ format
#############################################

##### Set pipeline parameters ###############
min_map_qual=2 # minimal first-pass mapping quality cutoff for slamdunk filtering
min_aln_identity=0.8 # minimal first-pass read alignment sequence identify cutoff for slamdunk filtering
min_edit_distance=-1 # minimal first-pass mapping edit distance cutoff for slamdunk filtering

min_cov_for_snv=10 # minimal variant coverage cutoff for Varscan calling
min_var_freq=0.2 # minimal variant frequency cutoff for Varscan calling

min_map_qual2=30 # minimal mapping quality cutoff for base-converted reads extraction
min_converted_base_count=2 # minimal converted base count cutoff for base-converted reads extraction

threads="24" # minimal CPU threads to use
debug="no" # Whether to keep intermediate files for debugging: "yes" or "no". Default: "no"
############################################

############################################
#Normally no need to change the code below #
############################################

adapters="$trimmomatic_dir/adapters/*.fa"
cat $adapters |sed 's/>/\n>/' > adapter.fa

mkdir tmp

# index reference sequence
if [[ ! -e ref.genome.fa.fai ]]
then
    $samtools_dir/samtools faidx ref.genome.fa
fi
if [[ ! -e ref.genome.dict ]]
then
    $java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary  -REFERENCE ref.genome.fa -OUTPUT ref.genome.dict
fi

# trim the reads
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33 $R1_reads $R2_reads $sample_id.R1.trimmed.PE.fq.gz $sample_id.R1.trimmed.SE.fq.gz $sample_id.R2.trimmed.PE.fq.gz $sample_id.R2.trimmed.SE.fq.gz ILLUMINACLIP:adapter.fa:2:30:10  SLIDINGWINDOW:5:20 MINLEN:36

# reverse complement R2 reads
$seqtk_dir/seqtk seq -r $sample_id.R2.trimmed.PE.fq.gz |gzip -c > $sample_id.R2.trimmed.revcom.PE.fq.gz

# first-pass mapping with ngm
$ngm_dir/ngm \
    -r $refseq \
    --qry1 $sample_id.R1.trimmed.PE.fq.gz \
    --qry2 $sample_id.R2.trimmed.revcom.PE.fq.gz \
    -t $threads \
    -n 1 \
    --strata \
    --bam \
    --slam-seq 2 \
    --no-progress \
    -o $sample_id.bam

$samtools_dir/samtools view -h -@ $threads $sample_id.bam | $samblaster_dir/samblaster | $samtools_dir/samtools sort -@ $threads -O bam -T $sample_id - >$sample_id.sort.bam
$samtools_dir/samtools index $sample_id.sort.bam

if [[ $debug == "no" ]]
then
    rm $sample_id.R1.trimmed.SE.fq.gz
    rm $sample_id.R2.trimmed.SE.fq.gz
    rm $sample_id.bam
fi

# slamdunk filtering
$slamdunk_dir/slamdunk filter -t $threads -mq $min_map_qual -mi $min_aln_identity -nm $min_edit_distance $sample_id.sort.bam -o  $sample_id.slamdunk_out

# bam2mpileup
$samtools_dir/samtools mpileup -B -A -f ref.genome.fa --output-QNAME ./$sample_id.slamdunk_out/$sample_id.sort_filtered.bam  > $sample_id.filtered.mpileup 

# Varscan SNV calling
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $varscan_dir/VarScan.jar mpileup2snp $sample_id.filtered.mpileup  --strand-filter 1 --output-vcf 1 --min-var-freq $min_var_freq --min-coverage $min_cov_for_snv -t $threads --variants 1 > $sample_id.vcf

$htslib_dir/bgzip $sample_id.vcf
$htslib_dir/tabix -p vcf $sample_id.vcf.gz

# vcf2reads
perl vcf2reads.pl -m $sample_id.filtered.mpileup -v $sample_id.vcf.gz -s $min_converted_base_count -q $min_map_qual2 -o $sample_id.read_list.txt 
$seqtk_dir/seqtk subseq $R1_reads $sample_id.read_list.txt |gzip -c  > $sample_id.base_converted.R1.fastq.gz
$seqtk_dir/seqtk subseq $R2_reads $sample_id.read_list.txt |gzip -c  > $sample_id.base_converted.R2.fastq.gz

# trim the reads
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33 $R1_reads $R2_reads $sample_id.base_converted.R1.trimmed.PE.fq.gz $sample_id.base_converted.R1.trimmed.SE.fq.gz $sample_id.base_converted.R2.trimmed.PE.fq.gz $sample_id.base_converted.R2.trimmed.SE.fq.gz ILLUMINACLIP:adapter.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:36

# reverse complement R2 reads
$seqtk_dir/seqtk seq -r $sample_id.base_converted.R2.trimmed.PE.fq.gz |gzip -c > $sample_id.base_converted.R2.trimmed.revcom.PE.fq.gz

# second-pass mapping for final visualization
$ngm_dir/ngm \
    -r $refseq \
    --qry1 $sample_id.base_converted.R1.trimmed.PE.fq.gz \
    --qry2 $sample_id.base_converted.R2.trimmed.revcom.PE.fq.gz \
    -t $threads \
    -n 1 \
    --strata \
    --bam \
    --slam-seq 2 \
    --no-progress \
    -o $sample_id.base_converted.bam

$samtools_dir/samtools index $sample_id.base_converted.bam
$samtools_dir/samtools view -h -@ $threads $sample_id.base_converted.bam | $samblaster_dir/samblaster | $samtools_dir/samtools sort -@ $threads -O bam -T $sample_id.base_converted - >$sample_id.base_converted.sort.bam
$samtools_dir/samtools index $sample_id.base_converted.sort.bam


if [[ $debug == "no" ]]
then
    rm $sample_id.base_converted.R1.trimmed.SE.fq.gz
    rm $sample_id.base_converted.R2.trimmed.SE.fq.gz
    rm -r $sample_id.sort.bam
    rm -r $sample_id.sort.bam.bai
    rm -r $sample_id.slamdunk_out
    rm $sample_id.vcf.gz
    rm $sample_id.vcf.gz.tbi
    rm $sample_id.filtered.mpileup
    rm $sample_id.base_converted.bam
    rm $sample_id.read_list.txt
fi

mkdir ${sample_id}_out
mv $sample_id.* ./${sample_id}_out
rm -r tmp
