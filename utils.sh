# add chr to chromosome name in a fasta file for nuclear chromosomes (tested on reference assembly from ensembl)
awk '{
	if ($0 ~ /^>([0-9]+|X|Y)( |:)/) {
        $0 = ">chr" substr($0, 2);
    }
    print
}' original.fasta > modified.fasta

# same as above but for annotation file (works for both gtf and gff3)
awk 'BEGIN {OFS="\t"} { if ($1 ~ /^([0-9]+|X|Y)$/ && substr($1, 1, 3) != "chr") $1 = "chr" $1; print }' original.gtf > modified.gtf

# get chr size from fasta index
cut -f1,2 input.fa.fai > input.size 

# extract first ~10k entries from a bam file. We are assuming header takes roughly 10 lines so use -n 10010. And output bam won't have exactly 10k entries.
samtools view -h original.bam | head -n 10010 | samtools view -b -o new.bam
sambamba index -t 8 new.bam
rm new_temp.bam

# extract reads that map to nuclear chromosomes from a bam file
# for mouse
samtools view -b input.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > output.bam
# for human
samtools view -b input.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > output.bam
