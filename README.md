# From Fastq to vcf

> Juliana Acosta-Uribe
Based on pipelines developed by [Khalid Mahmood](https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_gatk1/variant_calling_gatk1/) for Melbourne Bioinformatics, University of Melbourne (2021),  [Daniel Edward Deatherage](https://wikis.utexas.edu/display/bioiteam/Genome+Variant+Analysis+Course+2022) for the University of Texas (2022), [Mohammed Khalfan](https://learn.gencore.bio.nyu.edu/) for New York University, [Derek Caetano-Anolles](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) for the Broad institute (2023), [Griffith Lab](https://pmbio.org/module-04-germline/0004/02/02/Germline_SnvIndel_FilteringAnnotationReview/)



**Data:** \
Illumina HiSeq paired-end reads in FASTQ format from exomes. 
**Tools:** \
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://github.com/ewels/MultiQC), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [BWA-MEM](https://github.com/lh3/bwa), [Picard](https://broadinstitute.github.io/picard/), 

GATK4, Picard, Bcftools and jigv  
**Reference data:** \
GATK hg38 bundle of reference files downloaded from (ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/). 


## Tutorial contents table

* [Section 1: Perform Quality Control to fastq files](#section-1-perform-quality-control-to-fastq-files)
* [Section 2: Map raw mapped reads to reference genome](#section-2-map-raw-mapped-reads-to-reference-genome)
* [Section 3: Prepare analysis ready reads](#section-3-prepare-analysis-ready-reads)
* [Section 4: Variant calling](#section-4-variant-calling)
* [Section 5: Filter and prepare analysis ready variants](#section-5-filter-and-prepare-analysis-ready-variants)
* [Section 6: Exporting variant data and visualisation](#section-6-exporting-variant-data-and-visualisation)




## Section 1: Perform Quality Control to fastq files

### 1. Quality control with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://github.com/ewels/MultiQC)

The quality control of each .fastq will be analized using FastQC and then, we will use MultiQC to get an .html report of all the fastQC results.
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
for i in $SAMPLES
do 
java -jar trimmomatic-0.39.jar PE ./fastq/${i}_1.fastq.gz ./fastq/${i}_2.fastq.gz 
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

A. Create the BWA index files
```bash
bwa index hg38_bundle/Homo_sapiens_assembly38.fasta.gz
```
This will generate 5 additional files <reference.fasta.gz.amb>, <reference.fasta.gz.ann>, <reference.fasta.gz.bwt>, <reference.fasta.gz.pac>, <reference.fasta.gz.sa>\
Note: If the reference is greater than 2GB, you need to specify a different algorithm when building the BWA index, as follows: `bwa index -a bwtsw <reference.fasta>`

B. Align your fastq files
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

At the end of this step you should have UNSORTED `.bam` files

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
These duplicate reads are not informative and cannot be considered as evidence for or against a putative variant. For example, duplicates can arise during sample preparation e.g. library construction using PCR. Without this step, you risk having over-representation in your sequence of areas preferentially amplified during PCR. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates. \
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


###  3. Base quality recalibration
The last step of pre-processing mapped reads is the base quality score recalibration (BQSR) stage. The GATK tools detects systematic errors made by the sequencing machine while estimating the accuracy of each base. The systematic errors can be have various sources ranging from technical machine errors to the variability in the sequencing chemical reactions. The two step BQSR process applies machine learning to model the possible errors and adjust the base quality scores accordingly. Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly.
This allows us to get more accurate base qualities, which in turn improves the accuracy of our variant calls. More details [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-).

Get gatk running:
```bash
git clone https://github.com/broadinstitute/gatk.git
cd gatk/
./gradlew 
```
If you want to install gatk locally you can do `./gradlew localJar` [More info](https://github.com/broadinstitute/gatk)

**Step 1  - Build the model**
```bash
for i in 2 3
do
gatk --java-options "-Xmx7g" BaseRecalibrator \
-I ${i}.sorted.dup.bam \
-R reference/Homo_sapiens_assembly38.fasta \
--known-sites reference/dbsnp_146.hg38.vcf.gz \
--known-sites reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-O ${i}.recal_data.table
done
```
**#Step 2: Apply the model to adjust the base quality scores**
```bash
for i in $SAMPLES 
do
gatk ApplyBQSR \
-I ${i}.sorted.dup.bam \
-R reference/Homo_sapiens_assembly38.fasta \
--bqsr-recal-file ${i}.recal_data.table \
-O ${i}.sorted.dup.bqsr.bam 
done
```
We now have a pre-processed BAM file `sample.sorted.dup.bqsr.bam` ready for variant calling.

But before we proceed, let's take a detour and run some summary statistics of the alignment data and QC.

Get BQSR statistics
```bash
for i in $SAMPLES
do
gatk BaseRecalibrator \
-I ${i}.sorted.dup.bqsr.bam \
-R reference/Homo_sapiens_assembly38.fasta \
--known-sites reference/dbsnp_146.hg38.vcf.gz \
--known-sites reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-O ${i}.post-bqsr.recal_data.table

gatk AnalyzeCovariates \
-before ${i}.recal_data.table \
-after ${i}.post-bqsr.recal_data.table \
-plots ${i}.AnalyzeCovariates.pdf \

done
```
This will generate a .pdf file. More info [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037063232-AnalyzeCovariates)

Get BAM statistics and QC: \
The commands below uses FastQC and Picard to generate QC metrics followed by multiQC tools then aggregating the data to produce an HTML report.    
```bash
for i in $SAMPLES
do
java -jar picard.jar CollectMultipleMetrics \
R=reference/Homo_sapiens_assembly38.fasta \
I=${i}.sorted.dup.bqsr.bam \
O=${i}.sorted.dup.bqsr.CollectMultipleMetrics
done
multiqc . --filename post-bqsr
```
======================================May 26

## Section 3: Variant calling

The next step in the GATK best practices workflow is to proceed with the variant calling.

There are a couple of workflows to call variants using GATK4. Here we will follow the Genomic Variant Call Format (GVCF) workflow which is more suited for scalable variant calling i.e. allows incremental addition of samples for joint genotyping.

### 1. Apply HaplotypeCaller
HaplotypeCaller is the focal tool within GATK4 to simultaneously call germline SNVs and small Indels using local de-novo assembly of haplotype regions.

!!! Algorithm
    Briefly, the HaplotypeCaller works by:
    1. Identify active regions or regions with evidence of variations.
    2. Re-asssemble the active regions.
    3. Re-align active region reads to the assembled regions to identify allele.
    More details about the HaplotypeCaller can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller).

```bash
EXOME_TARGETS=file
for i in $SAMPLES
do
gatk HaplotypeCaller \
-I ${i}.sorted.dup.bqsr.bam \
-R reference/Homo_sapiens_assembly38.fasta \
-ERC GVCF \
#-L $EXOME_TARGETS \
-O ${i}.g.vcf.gz
done
```

The output of this step is a GVCF file. The format for the GVCF file is similar to a VCF file. The key difference is that the GVCF file contains records for each sequenced genomic coordinate. The `--emit-ref-confidence` or `-ERC` parameter lets you select a method to summarise confidence in the genomic site being homozygous-reference. The option `-ERC GVCF` is more efficient and recommended for large samples and therefore more scalable.

### 2. Apply CombineGVCFs
The CombineGVCFs tool is applied to combine multiple single sample GVCF files to merge these in to a single multi-sample GVCF file.

Create a text file containing the all GVCFs you want to combine:
```bash
ls *.vcf.gz > gvcfs.list
```
Merge the GVCFs into a single file
```bash
gatk CombineGVCFs \
-R reference/Homo_sapiens_assembly38.fasta \
-V gvcfs.list \
#-L $EXOME_TARGETS \
-O merged_gvcf.g.vcf.gz
```
Now that we have a merged GVCF file, we are ready to perform genotyping.

### 3. Apply GenotypeGVCFs
```bash
gatk GenotypeGVCFs \
-R reference/Homo_sapiens_assembly38.fasta \
-V merged_gvcf.g.vcf.gz
#-L $EXOME_TARGETS \
-O cohort.vcf.gz
```
??? Information
    An alternative to CombineGVCFs is [GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport), which is more efficient for sample numbers and stores the content in a a GenomicsDB data store. Therefore, CombineGVCFs could be slow and inefficient for more than a few samples. A possible work around is to split up the tasks per interval regions such as chromosomes.

## Section 4: Filter and prepare analysis ready variants

### 1. Variant Quality Score Recalibration

The raw VCF file from the previous step (`cohort.vcf.gz`) contains 10467 variants. Not all of these are real, therefore, the aim of this step is to filter out artifacts or false positive variants. The GATK Best Practices recommends filtering germline variant callsets with [VQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-). 
A second filtering strategy is called "Hard filtering", which is useful when the data cannot support VQSR or when an analysis requires manual filtering.  More detail [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering)

The Variant Quality Score Recalibration or the VQSR strategy is a two step process (1) the first step builds a model that describes how variant metric or quality measures co-vary with the known variants in the training set. (2) The second step then ranks each variant according to the target sensitivity cutoff and applies a filter expression.

**Split into SNPs and INDELs**
```bash
for i in SNP INDEL
do
gatk SelectVariants 
-R reference/Homo_sapiens_assembly38.fasta 
-V cohort.vcf.gz 
--select-type-to-include ${i} 
-O cohort.${i}.vcf.gz 
done
```
You should end with two new files `cohort.SNP.vcf.gz` and `cohort.INDEL.vcf.gz`

**Step 1: Build SNP model**
```bash
gatk VariantRecalibrator 
-R reference/Homo_sapiens_assembly38.fasta 
-V cohort.SNP.vcf.gz 
-O cohort.SNP.recal 
-tranche 100.0 
-tranche 99.9 
-tranche 99.0 
-tranche 90.0 
--tranches-file cohort.SNP.tranches 
--trust-all-polymorphic 
--mode SNP 
--max-gaussians 4 
--resource hapmap,known=false,training=true,truth=true,prior=15.0:reference/hapmap_3.3.hg38.vcf.gz 
--resource omni,known=false,training=true,truth=true,prior=12.0:reference/1000G_omni2.5.hg38.vcf.gz 
--resource 1000G,known=false,training=true,truth=false,prior=10.0:reference/1000G_phase1.snps.high_confidence.hg38.vcf.gz 
--resource dbsnp,known=true,training=false,truth=false,prior=2.0:reference/dbsnp_138.hg38.vcf.gz 
-an QD 
-an MQ 
-an MQRankSum 
-an ReadPosRankSum 
-an FS 
-an SOR 
--rscript-file recalibrate_SNP_plots.R
```
Note: These parameters are for exome data. 

**Step 2: Apply recalibration to SNPs**
```bash
gatk ApplyVQSR 
-R reference/Homo_sapiens_assembly38.fasta 
-V cohort.SNP.vcf.gz 
-O cohort.SNP.vqsr.vcf.gz 
--mode SNP 
--truth-sensitivity-filter-level 99.0 
--tranches-file cohort.SNP.tranches 
--recal-file cohort.SNP.recal 
```
**Step 3: Build Indel recalibration model**
```bash
gatk VariantRecalibrator 
-R reference/Homo_sapiens_assembly38.fasta \
-V cohort.INDEL.vcf.gz \
-O cohort.INDEL.recal \
-tranche 100.0 
-tranche 99.9 
-tranche 99.0 
-tranche 90.0
--tranches-file cohort.INDEL.tranches
--mode INDEL 
--max-gaussians 4 
--resource mills,known=false,training=true,truth=true,prior=12.0:reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 
--resource dbsnp,known=true,training=false,truth=false,prior=2.0:reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz 
-an QD 
-an FS 
-an SOR 
-an MQRankSum 
-an ReadPosRankSum 
--rscript-file recalibrate_INDEL_plots.R  
```
**Step 4: Apply recalibration to INDELs**
```bash
gatk ApplyVQSR 
-R reference/Homo_sapiens_assembly38.fasta \
-V cohort.INDEL.vcf.gz \
-O cohort.INDEL.vqsr.vcf.gz  \
--mode INDEL 
--truth-sensitivity-filter-level 99.0 
--tranches-file cohort.INDEL.tranches 
--recal-file cohort.INDEL.recal 
```
