#!/bin/sh
# *******************************************
# Script to perform TN seq variant calling
# using a matched paired Tumor+normal sample with fastq
# files named normal_1.fastq.gz, normal_2.fastq.gz
# tumor_1.fastq.gz, tumor_2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
fastq_folder="/home/pipeline/samples"
tumor_fastq_1="tumor_1.fastq.gz"
tumor_fastq_2="tumor_2.fastq.gz" #If using Illumina paired data
tumor_sample="tumor_sample_name"
tumor_group="tumor_read_group_name"
normal_fastq_1="normal_1.fastq.gz"
normal_fastq_2="normal_2.fastq.gz" #If using Illumina paired data
normal_sample="normal_sample_name"
normal_group="normal_read_group_name"
platform="ILLUMINA"

# Update with the location of the reference data files
fasta="/home/regression/references/b37/human_g1k_v37_decoy.fasta"
dbsnp="/home/regression/references/b37/dbsnp_138.b37.vcf.gz"
known_Mills_indels="/home/regression/references/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
known_1000G_indels="/home/regression/references/b37/1000G_phase1.indels.b37.vcf.gz"

# Update with the location of the Sentieon software package and license file
export SENTIEON_INSTALL_DIR=/home/release/sentieon-genomics-201911.01
export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic

# Other settings
nt=$(nproc) #number of threads to use in computation, set to number of cores in the server
workdir="$PWD/test/TNseq" #Determine where the output files will be stored

# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

# ******************************************
# 1a. Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$tumor_group\tSM:$tumor_sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$tumor_fastq_1 $fastq_folder/$tumor_fastq_2 || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o tumor_sorted.bam -t $nt --sam2bam -i -
# ******************************************
# 1b. Mapping reads with BWA-MEM, sorting for normal sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$normal_group\tSM:$normal_sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$normal_fastq_1 $fastq_folder/$normal_fastq_2 || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o normal_sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_sorted.bam --algo MeanQualityByCycle tumor_mq_metrics.txt --algo QualDistribution tumor_qd_metrics.txt --algo GCBias --summary tumor_gc_summary.txt tumor_gc_metrics.txt --algo AlignmentStat --adapter_seq '' tumor_aln_metrics.txt --algo InsertSizeMetricAlgo tumor_is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o tumor_gc-report.pdf tumor_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o tumor_qd-report.pdf tumor_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o tumor_mq-report.pdf tumor_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o tumor_is-report.pdf tumor_is_metrics.txt

# ******************************************
# 2b. Metrics for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i normal_sorted.bam --algo MeanQualityByCycle normal_mq_metrics.txt --algo QualDistribution normal_qd_metrics.txt --algo GCBias --summary normal_gc_summary.txt normal_gc_metrics.txt --algo AlignmentStat --adapter_seq '' normal_aln_metrics.txt --algo InsertSizeMetricAlgo normal_is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o normal_gc-report.pdf normal_gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o normal_qd-report.pdf normal_qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o normal_mq-report.pdf normal_mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o normal_is-report.pdf normal_is_metrics.txt

# ******************************************
# 3a. Remove Duplicate Reads for tumor
# sample. It is possible
# to mark instead of remove duplicates
# by ommiting the --rmdup option in Dedup
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i tumor_sorted.bam --algo LocusCollector --fun score_info tumor_score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i tumor_sorted.bam --algo Dedup --rmdup --score_info tumor_score.txt --metrics tumor_dedup_metrics.txt tumor_deduped.bam 
# ******************************************
# 3b. Remove Duplicate Reads for normal
# sample. It is possible
# to mark instead of remove duplicates
# by ommiting the --rmdup option in Dedup
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i normal_sorted.bam --algo LocusCollector --fun score_info normal_score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i normal_sorted.bam --algo Dedup --rmdup --score_info normal_score.txt --metrics normal_dedup_metrics.txt normal_deduped.bam 

# ******************************************
# 2a. Coverage metrics for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i normal_deduped.bam --algo CoverageMetrics normal_coverage_metrics
# ******************************************
# 2b. Coverage metrics for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_deduped.bam --algo CoverageMetrics tumor_coverage_metrics

# ******************************************
# 4a. Indel realigner for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels tumor_realigned.bam
# ******************************************
# 4b. Indel realigner for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i normal_deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels normal_realigned.bam

# ******************************************
# 5a. Base recalibration for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels tumor_recal_data.table
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -q tumor_recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels tumor_recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt --algo QualCal --plot --before tumor_recal_data.table --after tumor_recal_data.table.post tumor_recal.csv
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o tumor_recal_plots.pdf tumor_recal.csv
# ******************************************
# 5b. Base recalibration for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i normal_realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels normal_recal_data.table
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i normal_realigned.bam -q normal_recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels normal_recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt --algo QualCal --plot --before normal_recal_data.table --after normal_recal_data.table.post normal_recal.csv
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o normal_recal_plots.pdf normal_recal.csv

# ******************************************
# 6. Corealignment of tumor and normal
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -i normal_realigned.bam -q tumor_recal_data.table -q normal_recal_data.table --algo Realigner -k $known_Mills_indels -k $known_1000G_indels tn_corealigned.bam

# ******************************************
# 7. Somatic Variant Calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tn_corealigned.bam --algo TNsnv --tumor_sample $tumor_sample --normal_sample $normal_sample --dbsnp $dbsnp --call_stats_out output-call.stats output-tnsnv.vcf.gz
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i tn_corealigned.bam --algo TNhaplotyper --tumor_sample $tumor_sample --normal_sample $normal_sample --dbsnp $dbsnp output-tnhaplotyper.vcf.gz
