#!/bin/sh
# *******************************************
# Script to perform joint calling on 3 samples
# with fastq files named sample<i>_1.fastq.gz
# and sample<i>_2.fastq.gz
# *******************************************

# Update with the fullpath location and suffix of the fastq files
fastq_folder="/home/pipeline/samples"
fastq_1_suffix="1.fastq.gz"
fastq_2_suffix="2.fastq.gz" #If using Illumina paired data
platform="ILLUMINA" #platform

# Update with the location of the reference data files
regions="/home/regression/references/b37/TruSeq_exome_targeted_regions.b37.bed"
fasta="/home/regression/references/b37/human_g1k_v37_decoy.fasta"
dbsnp="/home/regression/references/b37/dbsnp_138.b37.vcf.gz"
known_Mills_indels="/home/regression/references/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
known_1000G_indels="/home/regression/references/b37/1000G_phase1.indels.b37.vcf.gz"

# Update with the location of the Sentieon software package and license file
export SENTIEON_INSTALL_DIR=/home/release/sentieon-genomics-201911.01
export SENTIEON_LICENSE=/home/Licenses/Sentieon.lic

# Other settings
nt=$(nproc) #number of threads to use in computation, set to number of cores in the server
workdir="$PWD/test/DNAseq_joint" #Determine where the output files will be stored

# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1

# ******************************************
# 0. Process all samples independently
# ******************************************
#update with the prefix of the fastq files
for sample in sample1 sample2 sample3; do
  group=$sample
  mkdir $workdir/$sample
  cd $workdir/$sample

  # ******************************************
  # 1. Mapping reads with BWA-MEM, sorting
  # ******************************************
  ( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/${sample}_$fastq_1_suffix $fastq_folder/${sample}_$fastq_2_suffix || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort -r $fasta -o sorted.bam -t $nt --sam2bam -i -
  
  # ******************************************
  # 2. Metrics
  # ******************************************
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
  $SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o gc-report.pdf gc_metrics.txt
  $SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qd-report.pdf qd_metrics.txt
  $SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o mq-report.pdf mq_metrics.txt
  $SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o is-report.pdf is_metrics.txt
 
  # ******************************************
  # 3. Remove Duplicate Reads. It is possible
  # to mark instead of remove duplicates
  # by ommiting the --rmdup option in Dedup
  # ******************************************
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i sorted.bam --algo LocusCollector --fun score_info score.txt
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam 
  
  # ******************************************
  # 2a. Coverage metrics
  # ******************************************
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --algo CoverageMetrics coverage_metrics

  # ******************************************
  # 5. Base recalibration
  # ******************************************
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
  $SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o recal_plots.pdf recal.csv
  
  # ******************************************
  # 6b. HC Variant caller for GVCF
  # ******************************************
  $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta --interval $regions -t $nt -i deduped.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=30 --call_conf=30 --emit_mode gvcf output-hc.g.vcf.gz
  gvcf_inputs="$gvcf_inputs -v $workdir/$sample/output-hc.g.vcf.gz"
done

# ******************************************
# Perform the joint calling 
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta --algo GVCFtyper $gvcf_inputs -d $dbsnp $workdir/output-jc.vcf.gz
