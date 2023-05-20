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


###  3. Base quality recalibration
The last step of pre-processing mapped reads is the base quality score recalibration (BQSR) stage. The GATK tools detects systematic errors made by the sequencing machine while estimating the accuracy of each base. The systematic errors can be have various sources ranging from technical machine errors to the variability in the sequencing chemical reactions. The two step BQSR process applies machine learning to model the possible errors and adjust the base quality scores accordingly. More details [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-).

```bash
# lets go to the home directory again
cd

# step 1  - Build the model
gatk --java-options "-Xmx7g" BaseRecalibrator \
    -I output/NA12878.sort.dup.bam \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    --known-sites reference/hg38/dbsnp_146.hg38.vcf.gz \
    -O output/recal_data.table
```

```bash
# step 2: Apply the model to adjust the base quality scores
gatk --java-options "-Xmx7g" ApplyBQSR \
    -I output/NA12878.sort.dup.bam \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    --bqsr-recal-file output/recal_data.table \
    -O output/NA12878.sort.dup.bqsr.bam
```

!!! note
    In a workflow such as this it is a good practice to give the output files an appropriate name. In this case, we are appending the workflow step details to the filenames. For example, append `dup` after running the mark duplicates step.

We now have a pre-processed BAM file (`#!bash NA12878.sort.dup.bqsr.bam`) ready for variant calling.

But before we proceed, let's take a detour and run some summary statistics of the alignment data and QC.

??? example "**BAM statistics and QC** "
    The commands below uses FastQC and Picard to generate QC metrics followed by multiQC tools then aggregating the data to produce an HTML report.    

    ```bash
    # FastQC
    fastqc data/NA12878.chr20.region_1.fastq.gz \
    data/NA12878.chr20.region_2.fastq.gz \
    -o output/

    # CollectInsertSizeMetrics
    picard CollectMultipleMetrics \
    R=reference/hg38/Homo_sapiens_assembly38.fasta \
    I=output/NA12878.sort.dup.bqsr.bam \
    O=output/NA12878.sort.dup.bqsr.CollectMultipleMetrics

    # MultiQC
    multiqc output/. -o output/.    
    ```
    We have precomputed this and the resulting MultiQC report is [here](files/multiqc_report.html){:target="_blank"}.


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
gatk --java-options "-Xmx7g" HaplotypeCaller \
    -I output/NA12878.sort.dup.bqsr.bam \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -ERC GVCF \
    -L chr20 \
    -O output/NA12878.g.vcf.gz
```

The output of this step is a GVCF file. The format for the GVCF file is similar to a VCF file. The key difference is that the GVCF file contains records for each sequenced genomic coordinate. The `#!bash --emit-ref-confidence` or `#!bash -ERC` parameter lets you select a method to summarise confidence in the genomic site being homozygous-reference. The option `#!bash -ERC GVCF` is more efficient and recommended for large samples and therefore more scalable.

### 2. Apply CombineGVCFs
The CombineGVCFs tool is applied to combine multiple single sample GVCF files to merge these in to a single multi-sample GVCF file.

We have pre-processed two additional samples (NA12891 and NA12892) up to the HaplotypeCaller step (above). Lets first copy the gvcf files to the output directory.

```bash
#lets make sure that we are in the apropriate directory
cd

cp /mnt/shared_data/NA12891.g.vcf.gz* output/.
cp /mnt/shared_data/NA12892.g.vcf.gz* output/.

```

```bash
gatk --java-options "-Xmx7g" CombineGVCFs \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -V output/NA12878.g.vcf.gz \
    -V output/NA12891.g.vcf.gz \
    -V output/NA12892.g.vcf.gz \
    -L chr20 \
    -O output/cohort.g.vcf.gz
```

??? example "Lets look at the combined gvcf file"
    ```bash
    less output/cohort.g.vcf.gz
    ```
    Work your way down to the variant records? How many samples do you see in the VCF file?
    Hint: look at the header row.

Now that we have a merged GVCF file, we are ready to perform genotyping.

### 3. Apply GenotypeGVCFs
GenotypeGVCFs

```bash
gatk --java-options "-Xmx7g" GenotypeGVCFs \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -V output/cohort.g.vcf.gz \
    -L chr20 \
    -O output/output.vcf.gz
```

??? Information
    An alternative to CombineGVCFs is [GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport), which is more efficient for sample numbers and stores the content in a a GenomicsDB data store. Therefore, CombineGVCFs could be slow and inefficient for more than a few samples. A possible work around is to split up the tasks per interval regions such as chromosomes.


??? example "Visualisations: VCF file"
    Screenshot from output.vcf.gz
    ![fig1](./media/fig1.png)


## Section 4: Filter and prepare analysis ready variants

### 1. Variant Quality Score Recalibration

The raw VCF file from the previous step (`#!bash output.vcf.gz`) contains 10467 variants. Not all of these are real, therefore, the aim of this step is to filter out artifacts or false positive variants. GATK has provided different workflows for variant filtering. Here we will walk through the Variant Quality Score Recalibration or the VQSR strategy. VQSR is a two step process (1) the first step builds a model that describes how variant metric or quality measures co-vary with the known variants in the training set. (2) The second step then ranks each variant according to the target sensitivity cutoff and applies a filter expression.

```bash
#Step 1 - VariantRecalibrator
gatk --java-options "-Xmx7g" VariantRecalibrator \
    -V output/output.vcf.gz \
    --trust-all-polymorphic \
    -mode SNP \
    --max-gaussians 6 \
    --resource:hapmap,known=false,training=true,truth=true,prior=15 reference/hg38/hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=true,prior=12 reference/hg38/1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10 reference/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=7 reference/hg38/dbsnp_138.hg38.vcf.gz \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -O output/cohort_snps.recal \
    --tranches-file output/cohort_snps.tranches

#Step 2 - ApplyVQSR
gatk --java-options "-Xmx7g" ApplyVQSR \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -V output/output.vcf.gz \
    -O output/output.vqsr.vcf \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file output/cohort_snps.tranches \
    --recal-file output/cohort_snps.recal \
    -mode SNP
```

!!! CountVariants
    There are number of ways to count the variants in a VCF file. A very straight forward way using the GATK4 tools is as follows:
    ```bash
    gatk CountVariants -V output/output.vqsr.vcf
    ```

    ```
    Tool returned:
    10467
    ```

There are several protocols for filtering VCF files. We have walked throught the VQSR strategy above and for other options please visit this [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering).

??? example "Filtering strategy for a single sample VCF file"
    Consider the following method to filter a single sample VCF file. Here we will go through the Convolutional Neural Net based protocol to annotate and filter the VCF file.

    This is a two step process:

    (i) CNNScoreVariants will annotate the variant with pre-computed single-sample derived model scores in the INFO field CNN_1D (the neural network performs convolutions over the reference sequence surrounding the variant and combines those features with a multilayer perceptron on the variant annotations).

    ```bash
    gatk --java-options "-Xmx7g" CNNScoreVariants  \
       -R reference/hg38/Homo_sapiens_assembly38.fasta \
       -V output/output.vcf.gz \
       -O output/output.cnns.vcf
    ```

    (ii) FilterVariantTranches takes as input the percent sensitivities (0-100) to known sites to apply the filter. Variants with scores higher than for e.g. 99th percentile of variants in the resources pass through the filter and will have `PASS` in the filter. Others will have a filter values like 'CNN_1D_INDEL_Tranche_99.40_100.00' or 'CNN_1D_SNP_Tranche_99.95_100.00'.

    ```bash
    gatk --java-options "-Xmx7g" FilterVariantTranches \
        -V output/output.cnns.vcf \
        --resource reference/hg38/hapmap_3.3.hg38.vcf.gz \
        --resource reference/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        --info-key CNN_1D \
        --snp-tranche 99.95 \
        --indel-tranche 99.4 \
        -O output/output.cnns.cnnfilter.vcf
    ```

!!! Hint
    BCFtools is a useful tool to manipulate, filter and query VCF files. More details from [BCFtools](https://samtools.github.io/bcftools/). BCFtools can be combined with linux command line tools as well to summarise data. For example, the command below can used extract and print the 'FILTER' column from the VCF file.

    ```bash
        bcftools query -f'%FILTER\n' output/output.vqsr.vcf.gz
    ```

### 2. Additional filtering
The VariantFiltration tools is designed for hard-filtering variant calls based on custom quality criteria such as sequencing depth, mapping quality etc. The two parameters are the filter-name and filter-expression. The parameter filter-name is the name of the filter to be used in the FILTER column if the expression in filter-expression is true. In the example below, if the sequencing depth at the variant site (VCF field DP) is less than 10, the FILTER field will be populated with the value 'Low_depth10'. Users can add multiple filter expression/name combinations.

```bash
gatk --java-options "-Xmx7g" VariantFiltration \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -V output/output.vqsr.vcf \
    -O output/output.vqsr.varfilter.vcf \
    --filter-name "Low_depth10" \
    --filter-expression "DP < 10"
```

!!! question "Question: How many variants have a low sequencing depth (DP<10) in the file output.vqsr.varfilter.vcf."

    ??? answer
        ```bash
            bcftools query -f'%FILTER\n' output/output.vqsr.varfilter.vcf | sort | uniq -c
        ```

        ```
          6 Low_depth10
          2 Low_depth10;VQSRTrancheSNP99.00to99.90
          9 Low_depth10;VQSRTrancheSNP99.90to100.00
          9064 PASS
          1278 VQSRTrancheSNP99.00to99.90
          108 VQSRTrancheSNP99.90to100.00
        ```

### 3. Final analysis ready VCF file

Given we have a filter annotated VCF files (), we can now create an analysis ready VCF file.


!!! question "Question: Create a VCF file called `output/output.vqsr.varfilter.pass.vcf.gz` that contains only PASS variants? The input VCF file is `output/output.vqsr.varfilter.vcf`." Hint: try using the Bcftools application."

    ??? answer
        Use the bcftools to filter PASS variants.
        ```bash
        bcftools view -f 'PASS,.' -O vcf -o output/output.vqsr.varfilter.pass.vcf output/output.vqsr.varfilter.vcf            
        ```

        We have now created an analysis ready version of the VCF file. It is also a good practice to compress and index the file.

        ```bash
        bgzip -c output/output.vqsr.varfilter.pass.vcf > output/output.vqsr.varfilter.pass.vcf.gz
        tabix -p vcf output/output.vqsr.varfilter.pass.vcf.gz
        ```


## Section 5: Exporting variant data and visualisation
VCF files, although in tabular format, are not user friendly. We will go through a couple of ways to share share and visualise variant data. This is important for downstream analysis as well as sharing data. First, we will convert the VCF file in to a TSV file (ready for Excel for example) in a manner where we extract data fields of interest.

### 1. VariantsToTable
This GATK4 tool extracts fields of interest from each record in a VCF file. [VariantsToTable](https://gatk.broadinstitute.org/hc/en-us/articles/360056968292-VariantsToTable) can extract field from both the INFO and FORMAT columns in the VCF file.

!!! note
    VariantsToTable, by default, only extracts PASS or . (no filtering applied) variants. Use the `#!bash --show-filtered` parameter to show all variants.

```bash
gatk VariantsToTable \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -V output/output.vqsr.varfilter.pass.vcf.gz \
    -F CHROM -F POS -F FILTER -F TYPE -GF AD -GF DP \
    --show-filtered \
    -O output/output.vqsr.varfilter.pass.tsv
```

### 2. HTML report
Another useful method for sharing data is an interactive HTML file. This is suited for sharing a smaller subset of variants along with sequencing data. Here we will go through a simple example using the [jigv](https://github.com/brentp/jigv) tool.

![fig2](./media/fig2.png)

We will start with creating a subset of variants to report.

```bash
bcftools view output/output.vqsr.varfilter.pass.vcf.gz \
chr20:3822018-3999324 | \
bgzip -c > output/subset.vcf.gz

tabix -p vcf output/subset.vcf.gz
```

Now, we will call the jigv tool command to generate the report.
```bash
jigv --sample NA12878 \
--sites output/subset.vcf.gz \
--fasta reference/hg38/Homo_sapiens_assembly38.fasta \
output/NA12878.sort.dup.bqsr.bam > output/NA12878.jigv.html
```

Here is an example [report](files/NA12878.html){:target="_blank"} we created earlier.

<!-- ```bash
jigv --sample NA12878 --ped family.ped --sites output/subset.vcf.gz --fasta reference/hg38/Homo_sapiens_assembly38.fasta output/*.sort.dup.bqsr.bam > output/NA12878.jigv.html
``` -->
