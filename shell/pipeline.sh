#!/bin/sh
#SBATCH -N 2
#SBATCH -p RM-shared
#SBATCH -t 5:00:00
#SBATCH --ntasks-per-node=128

#############
#Trim adapter
#############
##Rep1
#Round 1
cutadapt \
-f fastq \
--match-read-wildcards \
--times 1 \
-e 0.1 \
-O 1 \
--quality-cutoff 6 \
-m 18 \
-a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-g CTTCCGATCTACAAGTT \
-g CTTCCGATCTTGGTCCT \
-A AACTTGTAGATCGGA \
-A AGGACCAAGATCGGA \
-A ACTTGTAGATCGGAA \
-A GGACCAAGATCGGAA \
-A CTTGTAGATCGGAAG \
-A GACCAAGATCGGAAG \
-A TTGTAGATCGGAAGA \
-A ACCAAGATCGGAAGA \
-A TGTAGATCGGAAGAG \
-A CCAAGATCGGAAGAG \
-A GTAGATCGGAAGAGC \
-A CAAGATCGGAAGAGC \
-A TAGATCGGAAGAGCG \
-A AAGATCGGAAGAGCG \
-A AGATCGGAAGAGCGT \
-A GATCGGAAGAGCGTC \
-A ATCGGAAGAGCGTCG \
-A TCGGAAGAGCGTCGT \
-A CGGAAGAGCGTCGTG \
-A GGAAGAGCGTCGTGT \
-o rep1.IP.umi.r1.fqTr.fq \
-p rep1.IP.umi.r2.fqTr.fq \
rep1.IP.umi.r1.fq \
rep1.IP.umi.r2.fq

#fastqc -t 2 --extract -k 7 rep1.IP.umi.r1.fqTr.fq -o quality
#fastqc -t 2 --extract -k 7 rep1.IP.umi.r2.fqTr.fq -o quality

#Round 2
cutadapt \
-f fastq \
--match-read-wildcards \
--times 1 \
-e 0.1 \
-O 5 \
--quality-cutoff 6 \
-m 18 \
-A AACTTGTAGATCGGA \
-A AGGACCAAGATCGGA \
-A ACTTGTAGATCGGAA \
-A GGACCAAGATCGGAA \
-A CTTGTAGATCGGAAG \
-A GACCAAGATCGGAAG \
-A TTGTAGATCGGAAGA \
-A ACCAAGATCGGAAGA \
-A TGTAGATCGGAAGAG \
-A CCAAGATCGGAAGAG \
-A GTAGATCGGAAGAGC \
-A CAAGATCGGAAGAGC \
-A TAGATCGGAAGAGCG \
-A AAGATCGGAAGAGCG \
-A AGATCGGAAGAGCGT \
-A GATCGGAAGAGCGTC \
-A ATCGGAAGAGCGTCG \
-A TCGGAAGAGCGTCGT \
-A CGGAAGAGCGTCGTG \
-A GGAAGAGCGTCGTGT \
-o rep1.IP.umi.r1.fqTrTr.fq \
-p rep1.IP.umi.r2.fqTrTr.fq \
rep1.IP.umi.r1.fqTr.fq \
rep1.IP.umi.r2.fqTr.fq

#fastqc -t 2 --extract -k 7 rep1.IP.umi.r1.fqTrTr.fq -o quality
#fastqc -t 2 --extract -k 7 rep1.IP.umi.r2.fqTrTr.fq -o quality

##Rep2
#Round 1
cutadapt \
-f fastq \
--match-read-wildcards \
--times 1 \
-e 0.1 \
-O 1 \
--quality-cutoff 6 \
-m 18 \
-a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-g CTTCCGATCTACAAGTT \
-g CTTCCGATCTTGGTCCT \
-A AACTTGTAGATCGGA \
-A AGGACCAAGATCGGA \
-A ACTTGTAGATCGGAA \
-A GGACCAAGATCGGAA \
-A CTTGTAGATCGGAAG \
-A GACCAAGATCGGAAG \
-A TTGTAGATCGGAAGA \
-A ACCAAGATCGGAAGA \
-A TGTAGATCGGAAGAG \
-A CCAAGATCGGAAGAG \
-A GTAGATCGGAAGAGC \
-A CAAGATCGGAAGAGC \
-A TAGATCGGAAGAGCG \
-A AAGATCGGAAGAGCG \
-A AGATCGGAAGAGCGT \
-A GATCGGAAGAGCGTC \
-A ATCGGAAGAGCGTCG \
-A TCGGAAGAGCGTCGT \
-A CGGAAGAGCGTCGTG \
-A GGAAGAGCGTCGTGT \
-o rep2.IP.umi.r1.fqTr.fq \
-p rep2.IP.umi.r2.fqTr.fq \
rep2.IP.umi.r1.fq \
rep2.IP.umi.r2.fq

#Round 2
cutadapt \
-f fastq \
--match-read-wildcards \
--times 1 \
-e 0.1 \
-O 5 \
--quality-cutoff 6 \
-m 18 \
-A AACTTGTAGATCGGA \
-A AGGACCAAGATCGGA \
-A ACTTGTAGATCGGAA \
-A GGACCAAGATCGGAA \
-A CTTGTAGATCGGAAG \
-A GACCAAGATCGGAAG \
-A TTGTAGATCGGAAGA \
-A ACCAAGATCGGAAGA \
-A TGTAGATCGGAAGAG \
-A CCAAGATCGGAAGAG \
-A GTAGATCGGAAGAGC \
-A CAAGATCGGAAGAGC \
-A TAGATCGGAAGAGCG \
-A AAGATCGGAAGAGCG \
-A AGATCGGAAGAGCGT \
-A GATCGGAAGAGCGTC \
-A ATCGGAAGAGCGTCG \
-A TCGGAAGAGCGTCGT \
-A CGGAAGAGCGTCGTG \
-A GGAAGAGCGTCGTGT \
-o rep2.IP.umi.r1.fqTrTr.fq \
-p rep2.IP.umi.r2.fqTrTr.fq \
rep2.IP.umi.r1.fqTr.fq \
rep2.IP.umi.r2.fqTr.fq

###########################
#Remove repetitive elements
###########################
#Generate the repeat index
STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeSAindexNbases 9 \
--genomeDir homo_sapiens_repbase_v2 \
--genomeFastaFiles homo_sapiens_repbase_fixed_v2.fasta

#Rep1
#Run alignment
STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir homo_sapiens_repbase_v2 \
--genomeLoad NoSharedMemory \
--alignEndsType EndToEnd \
--outSAMunmapped Within \
--outFilterMultimapNmax 30 \
--outFilterMultimapScoreRange 1 \
--outFileNamePrefix rep1.IP.umi.r1.fqTrTr.sorted.STAR \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outBAMcompression 10 \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--outSAMattrRGline ID:foo \
--outSAMattributes All \
--outSAMmode Full \
--outStd Log \
--readFilesIn rep1.IP.umi.r1.fqTrTr.fq rep1.IP.umi.r2.fqTrTr.fq

#Rename files
mv rep1.IP.umi.r1.fqTrTr.sorted.STARAligned.out.bam rep1.IP.umi.r1.fq.repeat-mapped.bam
mv rep1.IP.umi.r1.fqTrTr.sorted.STARUnmapped.out.mate1 rep1.IP.umi.r1.fq.repeat-unmapped.fq
mv rep1.IP.umi.r1.fqTrTr.sorted.STARUnmapped.out.mate2 rep1.IP.umi.r2.fq.repeat-unmapped.fq

#Rep2
#Run alignment
STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir homo_sapiens_repbase_v2 \
--genomeLoad NoSharedMemory \
--alignEndsType EndToEnd \
--outSAMunmapped Within \
--outFilterMultimapNmax 30 \
--outFilterMultimapScoreRange 1 \
--outFileNamePrefix rep2.IP.umi.r1.fqTrTr.sorted.STAR \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outBAMcompression 10 \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--outSAMattrRGline ID:foo \
--outSAMattributes All \
--outSAMmode Full \
--outStd Log \
--readFilesIn rep2.IP.umi.r1.fqTrTr.fq rep2.IP.umi.r2.fqTrTr.fq

#Rename files
mv rep2.IP.umi.r1.fqTrTr.sorted.STARAligned.out.bam rep2.IP.umi.r1.fq.repeat-mapped.bam
mv rep2.IP.umi.r1.fqTrTr.sorted.STARUnmapped.out.mate1 rep2.IP.umi.r1.fq.repeat-unmapped.fq
mv rep2.IP.umi.r1.fqTrTr.sorted.STARUnmapped.out.mate2 rep2.IP.umi.r2.fq.repeat-unmapped.fq

##############################
#Map non-repeats to the genome
##############################
#Generate index 
STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir genome_STARindex \
--sjdbGTFfile ENCFF159KBI.gtf \
--genomeFastaFiles GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
--sjdbOverhang 99

##Rep1
#Alignment
STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir genome_STARindex \
--genomeLoad NoSharedMemory \
--readFilesIn \
rep1.IP.umi.r1.fq.repeat-unmapped.fq \
rep1.IP.umi.r2.fq.repeat-unmapped.fq \
--outSAMunmapped Within \
--outFilterMultimapNmax 1 \
--outFilterMultimapScoreRange 1 \
--outFileNamePrefix rep1.IP.umi.r1.fq.genome-mapped \
--outSAMattributes All \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--outSAMattrRGline ID:foo \
--outStd Log \
--alignEndsType EndToEnd \
--outBAMcompression 10 \
--outSAMmode Full

#Re-name BAM: rename genome-mapped outputs
mv rep1.IP.umi.r1.fq.genome-mappedAligned.out.bam rep1.IP.umi.r1.fq.genome-mapped.bam

##Rep2
#Alignment
STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir genome_STARindex \
--genomeLoad NoSharedMemory \
--readFilesIn \
rep2.IP.umi.r1.fq.repeat-unmapped.fq \
rep2.IP.umi.r2.fq.repeat-unmapped.fq \
--outSAMunmapped Within \
--outFilterMultimapNmax 1 \
--outFilterMultimapScoreRange 1 \
--outFileNamePrefix rep2.IP.umi.r1.fq.genome-mapped \
--outSAMattributes All \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--outSAMattrRGline ID:foo \
--outStd Log \
--alignEndsType EndToEnd \
--outBAMcompression 10 \
--outSAMmode Full

#Re-name BAM: rename genome-mapped outputs
mv rep2.IP.umi.r1.fq.genome-mappedAligned.out.bam rep2.IP.umi.r1.fq.genome-mapped.bam

###############
#Name sort BAM
###############
#Sort output from STAR by name to ensure read pairs are adjacent.
##Rep1
samtools \
sort \
-n \
-o rep1.IP.umi.r1.fq.genome-mappedSo.bam \
rep1.IP.umi.r1.fq.genome-mapped.bam

##Rep2
samtools \
sort \
-n \
-o rep2.IP.umi.r1.fq.genome-mappedSo.bam \
rep2.IP.umi.r1.fq.genome-mapped.bam

######################
#Remove PCR duplicates
######################
##Rep1
python barcodecollapsepe.py \
-o rep1.IP.umi.r1.fq.genome-mappedSo.rmDup.bam \
-m rep1.IP.umi.r1.fq.genome-mappedSo.rmDup.metrics \
-b rep1.IP.umi.r1.fq.genome-mappedSo.bam

##Rep2
python barcodecollapsepe.py \
-o rep2.IP.umi.r1.fq.genome-mappedSo.rmDup.bam \
-m rep2.IP.umi.r1.fq.genome-mappedSo.rmDup.metrics \
-b rep2.IP.umi.r1.fq.genome-mappedSo.bam

########################
#Sort resulting bam file
########################
##Rep1
samtools sort -o rep1.IP.umi.r1.fq.genome-mappedSo.rmDupSo.bam rep1.IP.umi.r1.fq.genome-mappedSo.rmDup.bam
##Rep2
samtools sort -o rep2.IP.umi.r1.fq.genome-mappedSo.rmDupSo.bam rep2.IP.umi.r1.fq.genome-mappedSo.rmDup.bam

#################
#Merge replicates
#################
samtools merge pe_clip.fq.genomemappedSo.rmDupSo.merged.bam rep1.IP.umi.r1.fq.genome-mappedSo.rmDupSo.bam rep2.IP.umi.r1.fq.genome-mappedSo.rmDupSo.bam