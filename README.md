# From Fastq to vcf

> Juliana Acosta-Uribe
Based on pipelines developed by [Khalid Mahmood](https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_gatk1/variant_calling_gatk1/) for Melbourne Bioinformatics, University of Melbourne (2021),  [Daniel Edward Deatherage](https://wikis.utexas.edu/display/bioiteam/Genome+Variant+Analysis+Course+2022) for the University of Texas (2022), [Mohammed Khalfan](https://learn.gencore.bio.nyu.edu/) for New York University



**Data:** Illumina HiSeq paired-end reads in FASTQ format from exomes. 
**Tools:**
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://github.com/ewels/MultiQC), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [BWA-MEM](https://github.com/lh3/bwa), [Picard](https://broadinstitute.github.io/picard/), 

GATK4, Picard, Bcftools and jigv  
**Reference data:** GATK hg38 bundle of reference files downloaded from (ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/). 


## Tutorial contents table

* [Section 1: Perform Quality Control to fastq files](#section-1-perform-quality-control-to-fastq-files)
* [Section 2: Map raw mapped reads to reference genome](#section-2-map-raw-mapped-reads-to-reference-genome)
* [Section 3: Prepare analysis ready reads](#section-3-prepare-analysis-ready-reads)
* [Section 4: Variant calling](#section-4-variant-calling)
* [Section 5: Filter and prepare analysis ready variants](#section-5-filter-and-prepare-analysis-ready-variants)
* [Section 6: Exporting variant data and visualisation](#section-6-exporting-variant-data-and-visualisation)




## Section 1: Perform Quality Control to fastq files

### 1. Quality control with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://github.com/ewels/MultiQC)

Analize fastq quality and get a report of all the fastQC results with MultiQC
```bash
SAMPLES="1 2 3"

for i in $SAMPLES
do 
fastqc ${i}_1.fastq.gz ${i}_2.fastq.gz -o .
done 

 multiqc . --filename original_fastq
```
Check the multiqc_report.html file that is generated

### 2. Trim and remove adapters with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) 

```bash
for i in $SAMPLES ; do \
java -jar trimmomatic-0.39.jar PE ./fastq/${i}_1.fastq.gz ./fastq/${i}_2.fastq.gz  \
${i}_1.paired-trimmed.fastq.gz ${i}_1.unpaired-trimmed.fastq.gz \
${i}_2.paired-trimmed.fastq.gz ${i}_2.unpaired-trimmed.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ; done
```
`TruSeq3-PE.fa` contains the information about the adapters that were used in your sequencing.
You must have the `TruSeq3-PE.fa` file in the directory you are running this command.

After running trimmomatic you will have 4 new files <file>_1.paired-trimmed.fastq.gz, <file>_1.unpaired-trimmed.fastq.gz, <file>_2.paired-trimmed.fastq.gz and  <file>_2.unpaired-trimmed.fastq.gz

### 3. Check Fastq files again
```bash
for i in $SAMPLES
do 
fastqc ${i}_1.paired-trimmed.fastq.gz ${i}_1.unpaired-trimmed.fastq.gz ${i}_2.paired-trimmed.fastq.gz ${i}_2.unpaired-trimmed.fastq.gz -o .
done 

multiqc . --filename pair_trimmed
```


## Section 2: Map raw mapped reads to reference genome

### 1. Preparation and data import

Prepare your Fastq files, and trim adapters if necessary. 
Download the *reference data*

### 2. Align genome with [BWA-MEM](https://github.com/lh3/bwa) 

Create the BWA index files
```bash
bwa index hg38_bundle/Homo_sapiens_assembly38.fasta.gz
```
This will generate 5 additional files <reference.fasta.gz.amb>, <reference.fasta.gz.ann>, <reference.fasta.gz.bwt>, <reference.fasta.gz.pac>, <reference.fasta.gz.sa>\
Note: If the reference is greater than 2GB, you need to specify a different algorithm when building the BWA index, as follows: `bwa index -a bwtsw <reference.fasta>`

Align your fastq
```bash
for i in $SAMPLES
do 
bwa mem -M -t 20 -R "@RG\tID:GNA_${i}\tSM:${i}\tPL:ILLUMINA" \
./hg38_bundle/Homo_sapiens_assembly38.fasta.gz \
./fastq/paired_trimmed/${i}_1.paired-trimmed.fastq.gz \
./fastq/paired_trimmed/${i}_2.paired-trimmed.fastq.gz | \
samtools view -b -h -o ${i}.bam 
done
```
There are two parts to the command here. The first part uses BWA to perform the alignment and the second part take the output from BWA and uses Samtools to convert the output to the BAM format.

**BWA flags:**\
`mem` Is used when the query sequences are longer than 70bp \
`-M` This flag tells bwa to consider split reads as secondary, required for GATK variant calling\
`-t` Number of threads (multi-threading mode)\
`-R` <readgroup_info> Provide the readgroup as a string. The read group information is key for downstream GATK functionality. The GATK will not work without a read group tag.\

**Samtools flags:**\
`-b, --bam` Output in the BAM format\
`-h, --with-header` Include the header in the output. \
`o FILE, --output` File name\

At the end of this step you should have UNSORTED .bam files

You can see that the `.bam` files have been properly annotated with the Read Groups `samtools view -H 1.bam | grep "RG"`

## Section 2: Prepare analysis ready reads

### 1. Sort SAM/BAM
The alignment file `<name>.bam` is not sorted. Before proceeding, we should sort the BAM file using the [Picard](https://broadinstitute.github.io/picard/) tools.

```bash
for i in $SAMPLES
do
java -jar picard.jar SortSam \
I=${i}.bam  \
O=${i}.sorted.bam  \
CREATE_INDEX=True \
SORT_ORDER=coordinate
done   
```

The above command will create a coordinate sorted BAM file and an index `<name.bai>` file.

!!! Alignment statistics
    Given we now have a sorted BAM file, we can now generate some useful statistics. To do so we can use the `samtools flagstat` command. More details are available [here](http://www.htslib.org/doc/samtools-flagstat.html).

```bash
samtools flagstat ${i}.sorted.bam
```

### 2. Mark duplicate reads
The aim of this step is to locate and tag duplicate reads in the BAM file.\
hese duplicate reads are not informative and cannot be considered as evidence for or against a putative variant. For example, duplicates can arise during sample preparation e.g. library construction using PCR. Without this step, you risk having over-representation in your sequence of areas preferentially amplified during PCR. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates. \
For more details go to [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-).

```bash
for i in $SAMPLES
do
java -jar picard.jar MarkDuplicates \
I=${i}.sorted.bam   \
O=${i}.sorted.dup.bam  \
M=marked_dup_metrics_${i}.txt
done
```
Note that this step does not remove the duplicate reads, but rather flags them as such in the read’s SAM record. Downstream GATK tools will ignore reads flagged as duplicates by default.
