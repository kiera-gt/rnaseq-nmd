#!/bin/bash
scriptname=$0   # $0 is the name of the program
seqtype=targeted #either whole or targeted
batch=""
trim=no
kmersize=35
readLength=151
vars=~/data/scripts/variableslist.sh
###
#help function
HELP () {
	echo -e "\n\tUSAGE: $scriptname -b BATCH [-s (whole|targeted)] [-t (yes|no)] [-v VARIABLES_FILE]\n"
	echo -e "\tMANDATORY OPTIONS:\n"
	echo -e "\t\t-b BATCH\t\tName of the sequencing run, usually in the form of mmmYYYY\n\n"
	echo -e "\tADDITIONAL OPTIONS:\n"
	echo -e "\t\t-v VARIABLES_FILE\tLocation of file containing variable information for analysis scripts"
	echo -e "\t\t\t\t\tDEFAULT: ~/data/scripts/variableslist.sh\n"	
	echo -e "\t\t-s SEQTYPE\t\tType of RNA sequenced"  
	echo -e "\t\t\t\t\tDEFAULT: targeted\n"						
	echo -e "\t\t-t TRIM\t\t\tWhether the FASTQ files to be aligned have been trimmed"
	echo -e "\t\t\t\t\tDEFAULT: no\n"
	echo -e "\t\t-k KMER_SIZE\t\tkmer size to use in GATK HaplotypeCaller. Must be an integer."
	echo -e "\t\t\t\t\tDEFAULT: 35\n"
	echo -e "\t\t-l READ_LENGTH\t\tRead Length. Must be an integer."
	echo -e "\t\t\t\t\tDEFAULT: 151\n"
	exit 0
}
###
#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  HELP
fi
###
#Get Options from command line
while getopts :k:t:l:b:s:v:h opt
do
	case $opt in
	t)trim=$OPTARG;;
	b)batch=$OPTARG;;
	s)seqtype=$OPTARG;;
	k)kmersize=$OPTARG;;
	l)readLength=$OPTARG;;
	v)vars=$OPTARG;;
    h) HELP; exit;;
    \?) echo "ERROR: Invalid option: -$OPTARG" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
    :) echo "ERROR: Option -$OPTARG requires an argument" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
    *) echo "ERROR: Unimplemented option: -$OPTARG" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
  esac
done
shift $(($OPTIND -1))
###
#check to make sure variables file exists
if [ ! -f $vars ] 
then
    echo -e "ERROR: File $vars DOES NOT exist" >&2
	echo -e "Edit this scripts default or use the -v option to point the script to your file that contains variables for analysis scripts" >&2
    exit 1
fi
###
source $vars
###
#check if seqtype is whole or targeted, set seqfile to appropriate directory name
seqfile=""
if [ $seqtype == 'targeted' ]
then 
	seqfile=targetedRNA
elif [ $seqtype == 'whole' ]
then
	seqfile=wholemRNA
else
	echo -e "ERROR: \"$seqtype\" is not a valid argument for [-s] option. Valid arguments are \"targeted\" or \"whole\""
	exit 1
fi
#make sure kmer size is a positive integer
if ! [[ "$kmersize" =~ ^[0-9]+$ ]]
    then
        echo -e "ERROR: \"$kmersize\" is not a valid argument for [-k] option. kmer size must be an integer."
		exit 1
fi

if ! [[ "$readLength" =~ ^[0-9]+$ ]]
    then
        echo -e "ERROR: \"$readLength\" is not a valid argument for [-l] option. Read Length must be an integer."
		exit 1
fi
#check whether or not to trim, set fastq_dir appropriately
fastq_type=""
if [ $trim == 'yes' ]
then 
	fastq_type=fastq_files/trimmed
	echo "Will align TRIMMED fastq files"
elif [ $trim == 'no' ]
then
	fastq_type=fastq_files
	echo "Will align UNTRIMMED fastq files"
else
	echo -e "ERROR: \"$trim\" is not a valid argument for [-t] option. Valid arguments are \"yes\" or \"no\""
	exit 1
fi

#ensure directory containing samples folder and samples.txt exists
if [ ! -d $basedir/$seqfile/$batch ] 
then
    echo -e "ERROR: Directory $basedir/$seqfile/$batch DOES NOT exist" >&2
	echo -e "Check that option arguments -b (& -s if used) is correct" >&2
    exit 1
fi

#ensure samples.txt exists
if [ ! -f $basedir/$seqfile/$batch/samples.txt ] 
then
    echo -e "ERROR: File $basedir/$seqfile/$batch/samples.txt DOES NOT exist" >&2
	echo -e "Create the standard samples.txt file for this batch" >&2
    exit 1
fi

star_2pass_run () { 
	local sample=$1
	local dir=$2
	local readgroup=$3
	mkdir $dir/star_2pass #STAR needs outFileNamePrefix to already exist before running
	local step=star_2pass
	local pbsfile=$dir/$step/$step.pbs
	fastq_dir=$dir/$fastq_type
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for star run
	echo "$STAR --runThreadN 16 --genomeDir $genomeDir --sjdbGTFfile $annotations --readFilesIn $fastq_dir/"$fastqname"_L001_R1_001.fastq.gz,$fastq_dir/"$fastqname"_L002_R1_001.fastq.gz,$fastq_dir/"$fastqname"_L003_R1_001.fastq.gz,$fastq_dir/"$fastqname"_L004_R1_001.fastq.gz $fastq_dir/"$fastqname"_L001_R2_001.fastq.gz,$fastq_dir/"$fastqname"_L002_R2_001.fastq.gz,$fastq_dir/"$fastqname"_L003_R2_001.fastq.gz,$fastq_dir/"$fastqname"_L004_R2_001.fastq.gz --readFilesCommand zcat --outFileNamePrefix $dir/star_2pass/ --sjdbOverhang 100 --outSAMtype BAM SortedByCoordinate --outSAMattrRGline ID:$readgroup.1 PL:illumina PU:$readgroup.1 LB:$seqtype SM:$sample , ID:$readgroup.2 PL:illumina PU:$readgroup.2 LB:$seqtype SM:$sample , ID:$readgroup.3 PL:illumina PU:$readgroup.3 LB:$seqtype SM:$sample , ID:$readgroup.4 PL:illumina PU:$readgroup.4 LB:$seqtype SM:$sample --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMultimapNmax 1 --outSJfilterCountUniqueMin 5 5 5 5 --outSJfilterCountTotalMin 5 5 5 5 --outFilterType BySJout --outSJfilterReads Unique --quantMode GeneCounts --twopassMode Basic" >> $pbsfile #add commands for star run
	echo -e "mv $dir/star_2pass/Aligned.sortedByCoord.out.bam $dir/star_2pass/$sample.bam" >> $pbsfile
	echo -e "mv $dir/star_2pass/Log.final.out $dir/star_2pass/$sample.log.final.out" >> $pbsfile
	echo -e "mv $dir/star_2pass/SJ.out.tab $dir/star_2pass/$sample.sj.out.tab" >> $pbsfile
	echo -e "samtools index $dir/star_2pass/$sample.bam" >> $pbsfile
	echo -e "qsub $dir/star_2pass/counts.pbs\n" >> $pbsfile
	echo -e "qsub $dir/star_2pass/qorts.pbs\n" >> $pbsfile
	echo -e "cd $dir/psi/" >> $pbsfile 
	echo -e "qsub $dir/psi/psi.pbs" >> $pbsfile
	echo -e "cd $dir/splicing/" >> $pbsfile 
	echo -e "qsub $dir/splicing/removedups.pbs" >> $pbsfile
	echo -e "cd $dir/variant_calling/" >> $pbsfile #move into variant_calling dir so pbs outfile goes to here
	echo -e "qsub $dir/variant_calling/markdups.pbs\n" >> $pbsfile #submit next pbs file in chain
}

psi () {
	local sample=$1
	local dir=$2
	local step=psi
	mkdir $dir/psi
	local pbsfile=$dir/psi/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "$psi_script -s $sample -d $dir -l $readLength" >> $pbsfile
	echo -e "sed -i \"1i #\t#\t$sample\t$sample\t$sample\" $basedir/$seqfile/$batch/samples/$sampledir/psi/$sample.exonic_parts.psi" >> $pbsfile
}

readcounts () {
	local sample=$1
	local dir=$2
	local step=counts
	mkdir $dir/genecounts
	local pbsfile=$dir/star_2pass/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_little" >> $pbsfile #start pbs file for step
	echo -e "$genecounts_script $sample $dir $geneinfo" >> $pbsfile
	echo -e "ontarget=\`awk '{sum += \$2} END {print sum}' $dir/genecounts/$sample.274genes.ncbi.txt\`; total=\`head -n 9 $dir/star_2pass/$sample.log.final.out | tail -n 1 | awk '{print \$6}'\`; awk -v tot=\"\$total\" -v ont=\"\$ontarget\" 'BEGIN{print \"Reads on target: \"ont/tot}'" >> $pbsfile
}

qorts () {
	local sample=$1
	local dir=$2
	local readgroup=$3
	local step=qorts
	mkdir $basedir/$seqfile/$batch/metrics/qorts/$sample.$readgroup.1
	mkdir $basedir/$seqfile/$batch/metrics/qorts/$sample.$readgroup.2
	mkdir $basedir/$seqfile/$batch/metrics/qorts/$sample.$readgroup.3
	mkdir $basedir/$seqfile/$batch/metrics/qorts/$sample.$readgroup.4
	local pbsfile=$dir/star_2pass/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $qorts QC --verbose --stranded --readGroup $readgroup.1 $dir/star_2pass/$sample.bam $annotations $basedir/$seqfile/$batch/metrics/qorts/$sample.$readgroup.1/" >> $pbsfile
	echo -e "java -jar $qorts QC --verbose --stranded --readGroup $readgroup.2 $dir/star_2pass/$sample.bam $annotations $basedir/$seqfile/$batch/metrics/qorts/$sample.$readgroup.2/" >> $pbsfile
	echo -e "java -jar $qorts QC --verbose --stranded --readGroup $readgroup.3 $dir/star_2pass/$sample.bam $annotations $basedir/$seqfile/$batch/metrics/qorts/$sample.$readgroup.3/" >> $pbsfile
	echo -e "java -jar $qorts QC --verbose --stranded --readGroup $readgroup.4 $dir/star_2pass/$sample.bam $annotations $basedir/$seqfile/$batch/metrics/qorts/$sample.$readgroup.4/" >> $pbsfile
}

remove_dups () {
	local sample=$1
	local dir=$2
	local step=removedups
	local pbsfile=$dir/splicing/$step.pbs
	mkdir $dir/splicing
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $picard MarkDuplicates I=$dir/star_2pass/$sample.bam O=$dir/splicing/$sample.nodups.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true M=$dir/splicing/removedups_metrics.txt" >> $pbsfile
	echo -e "samtools view -h -u -f 0x2 $dir/splicing/$sample.nodups.bam | samtools sort -n - | samtools fastq -1 $dir/splicing/$sample.nodups.R1.fastq -2 $dir/splicing/$sample.nodups.R2.fastq -" >> $pbsfile
	echo -e "gzip $dir/splicing/$sample.nodups.R1.fastq" >> $pbsfile
	echo -e "gzip $dir/splicing/$sample.nodups.R2.fastq" >> $pbsfile
	echo -e "cd $dir/splicing/" >> $pbsfile
	echo -e "qsub $dir/splicing/star_remap.pbs" >> $pbsfile
}

splicing_remap () {
	local sample=$1
	local dir=$2
	local step=star_remap
	local pbsfile=$dir/splicing/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for star run
	echo "$STAR --runThreadN 16 --genomeDir $genomeDir --sjdbGTFfile $annotations --readFilesIn $dir/splicing/$sample.nodups.R1.fastq.gz $dir/splicing/$sample.nodups.R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix $dir/splicing/ --sjdbOverhang 100 --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMultimapNmax 1 --outSJfilterCountUniqueMin 5 5 5 5 --outSJfilterCountTotalMin 5 5 5 5 --outFilterType BySJout --outSJfilterReads Unique --quantMode GeneCounts --twopassMode Basic" >> $pbsfile #add commands for star run
	echo -e "mv $dir/splicing/Aligned.sortedByCoord.out.bam $dir/splicing/$sample.remap.nodups.bam" >> $pbsfile
	echo -e "samtools index $dir/splicing/$sample.remap.nodups.bam" >> $pbsfile
	echo -e "mv $dir/splicing/Log.final.out $dir/splicing/$sample.nodups.log.final.out" >> $pbsfile
	echo -e "mv $dir/splicing/SJ.out.tab $dir/splicing/$sample.nodups.sj.out.tab" >> $pbsfile
	echo -e "cd $dir/splicing/" >> $pbsfile
	echo -e "qsub $dir/splicing/splice_counts.pbs" >> $pbsfile
	echo -e "qsub $dir/splicing/counts_nodups.pbs" >> $pbsfile
}

readcounts_nodups () {
	local sample=$1
	local dir=$2
	local step=counts_nodups
	local pbsfile=$dir/splicing/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_little" >> $pbsfile #start pbs file for step
	echo -e "$genecounts_nodupsscript $sample $dir $geneinfo" >> $pbsfile
	echo -e "ontarget=\`awk '{sum += \$2} END {print sum}' $dir/genecounts/$sample.nodups.274genes.ncbi.txt\`; total=\`head -n 9 $dir/splicing/$sample.nodups.log.final.out | tail -n 1 | awk '{print \$6}'\`; awk -v tot=\"\$total\" -v ont=\"\$ontarget\" 'BEGIN{print \"Reads on target: \"ont/tot}'" >> $pbsfile
}

norm_label_splices () {
	local sample=$1
	local dir=$2
	local step=splice_counts
	local pbsfile=$dir/splicing/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_little" >> $pbsfile #start pbs file for star run
	echo -e "$ind_splicing $dir $sample $sortedintervals $introns" >> $pbsfile
}

mark_duplicates () {
	local sample=$1
	local dir=$2
	local step=markdups
	mkdir $dir/variant_calling
	mkdir $dir/variant_calling/temp
	local pbsfile=$dir/variant_calling/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $picard MarkDuplicates I=$dir/star_2pass/$sample.bam O=$dir/variant_calling/temp/$sample.dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=true M=$dir/variant_calling/markdups_metrics.txt" >> $pbsfile
	echo -e "cd $dir/variant_calling/" >> $pbsfile #move into variant_calling dir so pbs outfile goes to here
	echo -e "qsub $dir/variant_calling/splitn.pbs\n" >> $pbsfile
}

split_reads () {
	local sample=$1
	local dir=$2
	local step=splitn
	local pbsfile=$dir/variant_calling/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $gatk -T SplitNCigarReads -R $wholeGenomeFasta -I $dir/variant_calling/temp/$sample.dedup.bam -o $dir/variant_calling/temp/$sample.dedup.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS" >> $pbsfile
	echo -e "cd $dir/variant_calling/" >> $pbsfile #move into variant_calling dir so pbs outfile goes to here
	echo -e "qsub $dir/variant_calling/baserecal.pbs\n" >> $pbsfile
}

base_recalibration () {
	local sample=$1
	local dir=$2
	local step=baserecal
	local pbsfile=$dir/variant_calling/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $gatk -T BaseRecalibrator -R $wholeGenomeFasta -I $dir/variant_calling/temp/$sample.dedup.split.bam -knownSites $phase1snps -knownSites $millsIndels -knownSites $dbsnp -o $dir/variant_calling/recal_data.table" >> $pbsfile
	echo -e "java -jar $gatk -T PrintReads -R $wholeGenomeFasta -I $dir/variant_calling/temp/$sample.dedup.split.bam -BQSR $dir/variant_calling/recal_data.table -o $dir/variant_calling/temp/$sample.dedup.split.recal.bam" >> $pbsfile
	echo -e "cd $dir/variant_calling/" >> $pbsfile #move into variant_calling dir so pbs outfile goes to here
	echo -e "qsub $dir/variant_calling/vcf.pbs\n" >> $pbsfile
}

vcf () {
	local sample=$1
	local dir=$2
	local step=vcf
	local pbsfile=$dir/variant_calling/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $gatk -T HaplotypeCaller -R $wholeGenomeFasta -I $dir/variant_calling/temp/$sample.dedup.split.recal.bam --dbsnp $dbsnp -kmerSize $kmersize -dontUseSoftClippedBases -stand_call_conf 20.0 -L $intervals -o $dir/variant_calling/$sample.nodups.raw.snps.indels.vcf" >> $pbsfile
	echo -e "java -jar $gatk -T VariantFiltration -R $wholeGenomeFasta -V $dir/variant_calling/$sample.nodups.raw.snps.indels.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -filterName DP -filter \"DP < 11.0\" -G_filterName GQ -G_filter \"GQ < 5.0\" -G_filterName HET -G_filter \"isHet == 1 && DP > 10\" -o $dir/variant_calling/$sample.nodups.filtered.vcf" >> $pbsfile
	echo -e "perl $annovar $dir/variant_calling/$sample.nodups.filtered.vcf $tools/annovar/humandb/ -buildver hg38 -out $dir/variant_calling/$sample.nodups.filtered.annotated.vcf -remove -protocol refGene,cytoBand,genomicSuperDups,avsnp150,dbnsfp33a,clinvar_20170905,exac03,gnomad_exome -operation g,r,r,f,f,f,f,f -argument '-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs' -nastring . -vcfinput" >> $pbsfile
}

mkdir $basedir/$seqfile/$batch/metrics/qorts
while read sample sampledir readgroup fastqname
do
	dir=$basedir/$seqfile/$batch/samples/$sampledir
	
	star_2pass_run $sample $dir $readgroup
	remove_dups $sample $dir
	splicing_remap $sample $dir
	readcounts_nodups $sample $dir
	norm_label_splices $sample $dir
	psi $sample $dir
	readcounts $sample $dir
	qorts $sample $dir $readgroup
	mark_duplicates $sample $dir
	split_reads $sample $dir
	base_recalibration $sample $dir
	vcf $sample $dir
    
	cd $dir/star_2pass/
	qsub $dir/star_2pass/star_2pass.pbs
done <$basedir/$seqfile/$batch/samples.txt
