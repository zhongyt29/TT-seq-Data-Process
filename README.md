# TT-seq-Data-Process
TT-seq data process and downstream analysis

## 1.Map
STAR比对
去除低质量reads，保留Proper_paired reads
去重

```bash
THREADS=$1
HUMANIDX=$2
FQ1=$3
FQ2=$4
SAMPLE=$5
ALIGNDIR=$6
TMPDIR=$7
SPIKEDIR=$8
SPIKEIDX=$9

#map TT-seq reads
STAR --runThreadN ${THREADS} --runMode alignReads --genomeDir ${HUMANIDX} --readFilesIn ${FQ1} ${FQ2} --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMunmapped None --outSAMattrRGline ID:${SAMPLE} PU:${SAMPLE} SM:${SAMPLE} LB:unknown PL:illumina --outSAMtype BAM Unsorted --outTmpDir ${TMPDIR}${SAMPLE} --outFileNamePrefix ${ALIGNDIR}${SAMPLE}
samtools view -q 10 -f 0x2 -b ${ALIGNDIR}${SAMPLE}Aligned.out.bam > ${ALIGNDIR}${SAMPLE}.pe.q10.bam
samtools sort --threads ${THREADS} -o ${ALIGNDIR}${SAMPLE}.sorted.bam ${ALIGNDIR}${SAMPLE}.pe.q10.bam
picard MarkDuplicates INPUT=${ALIGNDIR}${SAMPLE}.sorted.bam OUTPUT=${ALIGNDIR}${SAMPLE}.sorted.rmdup.bam METRICS_FILE=${ALIGNDIR}${SAMPLE}.sorted.marked.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT TMP_DIR=${TMPDIR}${SAMPLE}
samtools index ${ALIGNDIR}${SAMPLE}.sorted.rmdup.bam
rm ${ALIGNDIR}${SAMPLE}.sorted.bam

#map Spike-in reads
STAR --runThreadN ${THREADS} --runMode alignReads --genomeDir ${SPIKEIDX} --readFilesIn ${FQ1} ${FQ2} --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMunmapped None --outSAMattrRGline ID:${SAMPLE} PU:${SAMPLE} SM:${SAMPLE} LB:unknown PL:illumina --outSAMtype BAM Unsorted --outTmpDir ${TMPDIR}${SAMPLE}.spike --outFileNamePrefix ${ALIGNDIR}${SAMPLE}
samtools view -q 10 -f 0x2 -b ${SPIKEDIR}${SAMPLE}Aligned.out.bam > ${SPIKEDIR}${SAMPLE}.pe.q10.bam
samtools sort --threads ${THREADS} -o ${SPIKEDIR}${SAMPLE}.sorted.bam ${SPIKEDIR}${SAMPLE}.pe.q10.bam
picard MarkDuplicates INPUT=${SPIKEDIR}${SAMPLE}.sorted.bam OUTPUT=${SPIKEDIR}${SAMPLE}.sorted.rmdup.bam METRICS_FILE=${SPIKEDIR}${SAMPLE}.sorted.marked.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT TMP_DIR=${TMPDIR}${SAMPLE}
samtools index ${SPIKEDIR}${SAMPLE}.sorted.rmdup.bam 
rm ${SPIKEDIR}${SAMPLE}.sorted.bam
```
## 2.calculate scale factors
根据spike-in计算scale factors
```
#spike-in reads number
samtools flagstat ${SPIKEDIR}${SAMPLE}.sorted.rmdup.bam
```
## 3.normalize TT-seq and make strand-specific bigwig files
--
skipNonCoveredRegions --binSize 1 --scaleFactor $scalefactor --normalizeUsing
None --filterRNAstrand forward

：--skipNonCoveredRegions
--binSize 1 --scaleFactor -$scalefactor --normalizeUsing None --filterRNAstrand
reverse

```bash
BIGWIGDIR=$1
TMPDIR=$2
ALIGNDIR=$3
SAMPLE=$4
THREADS=$5
SCALEFACTOR=$6

mkdir -p ${BIGWIGDIR}
mkdir -p ${TMPDIR}

BAM="${ALIGNDIR}${SAMPLE}.sorted.rmdup.bam"

#tmp bam files
BAMFOR="${TMPDIR}${SAMPLE}.fwd.bam"     # BAM file representing reads mapping to forward strand
BAMREV="${TMPDIR}${SAMPLE}.rev.bam"     # BAM file representing reads mapping to reverse strand
BAMFOR1="${TMPDIR}${SAMPLE}.fwd1.bam"
BAMFOR2="${TMPDIR}${SAMPLE}.fwd2.bam"
BAMREV1="${TMPDIR}${SAMPLE}.rev1.bam"
BAMREV2="${TMPDIR}${SAMPLE}.rev2.bam"

#result bigwig files
BIGWIG="${BIGWIGDIR}${SAMPLE}.bigwig"           # BIGWIG file representing all reads
BIGWIGFOR="${BIGWIGDIR}${SAMPLE}.for.bigwig"    # BIGWIG file representing reads mapping to forward strand
BIGWIGREV="${BIGWIGDIR}${SAMPLE}.rev.bigwig"    # BIGWIG file representing reads mapping to reverse strand

#create bigwig file for all reads
bamCoverage --scaleFactor $SCALEFACTOR -p $THREADS -b $BAM -o $BIGWIG

#create bigwig file for the forward strand
samtools view -b -f 128 -F 16 --threads $THREADS $BAM > $BAMFOR1
samtools view -b -f 80  --threads $THREADS $BAM > $BAMFOR2
samtools merge --threads $THREADS -f $BAMFOR $BAMFOR1 $BAMFOR2
samtools index $BAMFOR
bamCoverage --scaleFactor $SCALEFACTOR -p $THREADS -b $BAMFOR -o $BIGWIGFOR

#create bigwig file for the reverse strand
samtools view -b -f 144 --threads $THREADS $BAM > $BAMREV1
samtools view -b -f 64 -F 16 --threads $THREADS $BAM > $BAMREV2
samtools merge --threads $THREADS -f $BAMREV $BAMREV1 $BAMREV2
samtools index $BAMREV
bamCoverage --scaleFactor $SCALEFACTOR -p $THREADS -b $BAMREV -o $BIGWIGREV

#remove temporary files
rm $BAMFOR $BAMFFOR1 $BAMFFOR2 $BAMREV $BAMREV1 $BAMREV2

```



