#!/bin/bash

basedir=~/scratch/nmd_samples
tools=~/data/tools
scripts=~/data/scripts
humangenome=~/scratch/Homo_sapiens

genomeDir=~/scratch/genomeDir/illumina_grch38
wholeGenomeFasta=$humangenome/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
annotations=$humangenome/NCBI/GRCh38/Annotation/Genes/genes.gtf
exonparts=$humangenome/NCBI/GRCh38/Annotation/Genes/ncbi274.exonic_parts.forPSI.gff
sortedintervals=~/data/target_coordinates/nmd274genes_intervals_vsort.txt
geneinfo=~/data/target_coordinates/nmd274genes_gtf106_ensg.txt
introns=~/data/target_coordinates/condensed_introns_nooverlap.txt

STAR=$tools/star/bin/Linux_x86_64/STAR
qorts=$tools/qorts/QoRTs-STABLE.jar
picard=$tools/picard/picard.jar
bedtools=$tools/bedtools2/bin

gatk=$tools/gatk/GenomeAnalysisTK.jar
phase1snps=$humangenome/gatk_resource_bundle/1000G_phase1.snps.high_confidence.hg38.vcf
millsIndels=$humangenome/gatk_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf
dbsnp=$humangenome/gatk_resource_bundle/dbsnp_146.hg38.vcf
annovar=$tools/annovar/table_annovar.pl

genecounts_script=$scripts/genecounts_nmdgenes.pl
genecounts_nodupsscript=$scripts/genecounts.nodups.nmdgenes.pl
psi_script=$scripts/calc_psi.sh
ind_splicing=$scripts/nmdgenes_sjs.pl

pbshead_big="#PBS -l nodes=1:ppn=12\n#PBS -l mem=40gb\n#PBS -l walltime=12:00:00\n#PBS -q iw-shared-6\n#PBS -j oe\n\ncd \$PBS_O_WORKDIR"
pbshead_little="#PBS -l nodes=1:ppn=4\n#PBS -l mem=8gb\n#PBS -l walltime=12:00:00\n#PBS -q iw-shared-6\n#PBS -j oe\n\ncd \$PBS_O_WORKDIR"
