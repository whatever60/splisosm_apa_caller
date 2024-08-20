#!/bin/bash

samples=(s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12)

# input bam and output data
working_dir=/mnt/c/aws_data/data/arora_nc_2022

script_dir=/mnt/c/aws_data/SpatialSplicing/repeat_polyapipe/script

# external data
annot_gff3=/mnt/c/aws_data/data/ensembl/pub/release-98/gff3/mus_musculus/Mus_musculus.GRCm38.98.chr_prefix.gff3.gz
tissue_positions=/mnt/c/aws_data/data/10x/visium_fresh_frozen_adult_mouse_brain/spatial/tissue_positions_list.csv
raw_h5=/mnt/c/aws_data/data/10x/visium_fresh_frozen_adult_mouse_brain/Visium_Fresh_Frozen_Adult_Mouse_Brain_raw_feature_bc_matrix.h5

mkdir -p $working_dir/bam
mkdir -p $working_dir/macs3_out
mkdir -p $working_dir/umi_tools_out

for sample in ${samples[@]}
do
    (
        # STEP 0: split bam to positive strand and negative strand
        samtools view -bh -F 0x10 $working_dir/bam/$sample\_bam.bam > $working_dir/bam/$sample\_bam.plus.bam
        samtools view -bh -f 0x10 $working_dir/bam/$sample\_bam.bam > $working_dir/bam/$sample\_bam.minus.bam

        # STEP 0.5: index the new bam files (optional)
        sambamba index -t 8 $working_dir/bam/$sample\_bam.plus.bam
        sambamba index -t 8 $working_dir/bam/$sample\_bam.minus.bam
        
        # STEP 1: MACS3 call peaks
        # set gsize to hs for human and mm for mouse
        macs3 callpeak \
            -t $working_dir/bam/$sample\_bam.plus.bam \
            --format BAM \
            --gsize mm \
            --outdir $working_dir/macs3_out/ \
            --name $sample\_plus \
            --nomodel \
            --extsize 90 \
            --shift 0 \
            --call-summit \
            --nolambda
        macs3 callpeak \
            -t $working_dir/bam/$sample\_bam.minus.bam \
            --format BAM \
            --gsize mm \
            --outdir $working_dir/macs3_out/ \
            --name $sample\_minus \
            --nomodel \
            --extsize 90 \
            --shift 0 \
            --call-summit \
            --nolambda

        # STEP 1.5: Combine narrowpeak output and specify strand
        awk \
            'BEGIN {OFS="\t"} FNR==NR {print $1, $2, $3, $4, $5, "+", $7, $8, $9, $10; next} {print $1, $2, $3, $4, $5, "-", $7, $8, $9, $10}' \
            $working_dir/macs3_out/$sample\_plus_peaks.narrowPeak $working_dir/macs3_out/$sample\_minus_peaks.narrowPeak \
            > $working_dir/macs3_out/$sample\_peaks.narrowPeak

        # macs3 bdgcallpeak \
        #     -i $working_dir/bam/$sample\_bam.bam \
        #     -l 30 \
        #     -o $working_dir/macs3_out/$sample\_peaks.narrowPeak

        # STEP 2: Extract narrowpeak on genes 3' UTRs to another narrowPeak (for visual 
        # inspection) and a tsv (for downstream analysis)
        python $script_dir/utils.py \
            -i $working_dir/macs3_out/$sample\_peaks.narrowPeak \
            -a $annot_gff3 \
            -o $working_dir/macs3_out/$sample\_peaks.tsv \
            -n $working_dir/macs3_out/$sample\_peaks.utr.narrowPeak

        # STEP 3: Add tag of peak to bam file using the tsv we got in STEP 2, getting a new bam file
        python $script_dir/add_peak_tag.py \
            -p $working_dir/macs3_out/$sample\_peaks.tsv \
            -i $working_dir/bam/$sample\_bam.bam \
            -o $working_dir/bam/$sample\_bam.peaks.bam \
            --force

        # STEP 4: UMI_tools deduplication and count table generation
        # Our new tag added in STEP 3 is "UP", so we put that here. Also, I changed UM 
        # tag from UR (UMI reported by sequencer) to UB (UMI corrected by space ranger) on 7.29
        umi_tools count \
            --per-gene \
            --per-cell \
            --gene-tag UP \
            --cell-tag CB \
            --umi-tag UB \
            --extract-umi-method tag \
            -I $working_dir/bam/$sample\_bam.peaks.bam \
            -S $working_dir/umi_tools_out/$sample\_counts.tsv.gz

        # STEP 5: Extract relevant data for IsoSDE into one folder for each sample
        python $script_dir/construct_apa.py \
            -a $annot_gff3 \
            -f ensembl-gff3 \
            -p $working_dir/macs3_out/$sample\_peaks.tsv \
            -c $working_dir/umi_tools_out/$sample\_counts.tsv.gz \
            -s $tissue_positions \
            -d $raw_h5 \
            -o $working_dir/apa
        
    ) &
done

# Wait for all background processes to finish
wait
