#!/bin/bash

rootpath=/share/home/Z-RNA
Datapath=/share/home/Z-RNA/rawdata

cd ${rootpath}
mkdir fastqc fastp star plasmid

# Quality testing
cd ${Datapath}
ls -d *[0-9]|while read id;
do
	fastqc -t 8 -o ${rootpath}/fastqc/ ${Datapath}/${id}/${id}_1.fastq.gz
  fastqc -t 8 -o ${rootpath}/fastqc/ ${Datapath}/${id}/${id}_2.fastq.gz;
done

# remove adapter
ls -d *[0-9]|while read id;
do 
  fastp -i ${Datapath}/${id}/${id}_1.fastq.gz -o ${rootpath}/fastp/${id}_1.clean.fastq.gz -I ${Datapath}/${id}/${id}_2.fastq.gz -O ${rootpath}/fastp/${id}_2.clean.fastq.gz -q 25 -u 10  -l 50 -w 10 -p;
done

# STAR
# generate index
cd ${rootpath}/plasmid/
less ~/ref/human/gencode.v42.primary_assembly.annotation.gtf >>plasmid.gtf
less ~/ref/human/GRCh38.primary_assembly.genome.fa >>plasmid.fa

STAR --runMode genomeGenerate \
--runThreadN 10 \
--genomeDir ${rootpath}/plasmid/plasmid_index \
--genomeFastaFiles ${rootpath}/plasmid/plasmid.fa     \
--sjdbGTFfile ${rootpath}/plasmid/plasmid.gtf \
--sjdbOverhang 100

# mapping
cd ${rootpath}/fastp/
ls *_1.clean.fastq.gz|while read id; do
STAR --runThreadN 20 \
--genomeDir ${rootpath}/plasmid/plasmid_index \
--genomeLoad NoSharedMemory \
--outReadsUnmapped Fastx \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes All \
--outFilterType BySJout \
--outFilterMultimapNmax 1 \
--alignSJoverhangMin 1\
--alignSJDBoverhangMin 8 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.06 \ ##adjust according to editing characteristic
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--readFilesCommand zcat \
--readFilesIn ${id} ${id%_*}_2.fastq.gz \
--outFileNamePrefix ${rootpath}/star/${id%_*}_
done