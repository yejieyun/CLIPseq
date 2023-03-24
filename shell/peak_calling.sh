############
#Call peaks
############
#active environt
conda deactivate
conda activate peakachu

#index data 
samtools index ctrl.bam
samtools index pe_clip.fq.genomemappedSo.rmDupSo.merged.bam

#run peakachu
peakachu adaptive \
-M 200 \
-m 0.0 \
-f 2.0 \
-Q 0.05 \
-c ctrl.bam \
-t pe_clip.fq.genomemappedSo.rmDupSo.merged.bam

#convert peakachu peaks to homer peak file
cp initial_peaks.csv peaks.txt 
awk -F'\t' -v OFS='\t' 'NR ==1 {print "ID", $0; next} {print (NR-1), $0} ' peaks.txt > peaks_homer_1.txt
cut -f6,7,8,9 --complement peaks_homer_1.txt > homer_peaks.txt

#############
#Find motifs 
#############
#download hg38 genome and load module
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

#converst peak file to bed file
awk 'NR>1{
    if ($2 < $3) {
        print $1, $2, $3, "clip_peak_"NR-1,$8,$4;
        }
    else {
        print $1, $3, $2, "clip_peak_"NR-1,$8,$4;
        }
    }' initial_peaks.csv | tr ' ' '\t' > initial_peaks.bed

#extract genomic DNA of peaks
bedtools getfasta -fi hg38.fa -bed initial_peaks.bed > initial_peaks_seq.fa

#find motifs 
findMotifs.pl initial_peaks_seq.fa fasta homer_output/ 

###############
#Annotate peaks
###############
#download genome .gtf file
wget -O hg38_UCSC.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz
gunzip hg38_UCSC.gtf.gz

#annotate peaks
annotatePeaks.pl homer_peaks.txt hg38.fa -gtf hg38_UCSC.gtf > homer_peaks_annot.txt

#####################
#Functional Analysis
####################
#install and run RCAS
conda install bioconductor-rcas -c bioconda
Rscript RCAS.R

