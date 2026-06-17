# **Practical Assignment - Genomics**


## Section 1: Clinical Genomics - **15 marks**

1. List 5 annotations commonly added to a vcf that could be used to help filter down variants to a handful of diagnostic variants. Give a brief explanation of how each could help. **[2 marks]**  

2. You are helping a medical scientist find the disease-causing variant for a patient of Indigenous Australian descent. You find a variant in a candidate gene that matches the clinical phenotype, and look at its annotations to find its Gnomad frequency is 0.6%. What can the frequency of this variant tell you about the pathogenicity of the variant? **[1 mark]**

3. Referring to the case in question 2, you partner with ANU’s National Centre for Indigenous Genomics researchers (NCIG), who provide some more information on the variant in the disease-causing gene. They give you access to a database of 50,000 Indigenous individuals, and you count up the frequency of this variant and find it is now at a frequency of 8% in this new dataset. What can you infer from this? Choose one answer only. **[1 mark]**
	- a. Gnomad database has many more participants than 50,000, therefore the frequency is lower in Gnomad database than NCIG. 
	- b. Gnomad database does not have good representation of the Australian indigenous populations, hence the low frequency. Given this variant’s high frequency in the NCIG data, there is a good chance the variant belongs to benign private variation, rather than is pathogenic.
	- c. Regardless of its % in the NCIG database, the variant should be recommended to be classed as possibly pathogenic because it is <1% in Gnomad, since Gnomad is the gold standard resource. 
	- d. It is not possible to infer anything from this, since 50,000 is not a big enough sample size to calculate population frequency. 

4. What is a de novo mutation? **[1 mark]** 

5. What are the three hallmarks of a cnv event? Use either a deletion or duplication event as your example. **[1.5 marks]**


6. Is it possible to have an autosomal dominant homozygous mutation? Explain how. If you want to choose a few routine filters, would you set your variant curation software to filter by autosomal dominant homozygous mutation, or by autosomal dominant heterozygous? **[2.5 marks]**


7. You annotate 3 variants according to their Ensembl sequence ontology and want to perform some filtering based on effect on transcript. You find one is a “start retained” variant, one is “stop lost”, and the other is an “inframe deletion”. If you set your filtering to only keep variants with impact marked HIGH, which of them would be kept? **[1 mark]**

8. Refer to question 7 - based on effect on gene function, give the reasoning why you think each impact rating was attributed to each of the 3 ontology terms here. **[1.5 marks]**


9. When determing QC thresholds for determining if a sample should pass or fail, why is %bp covered at {threshold}X a better method than just overall mean depth? **[1 mark]**

10. Sex check, contamination check, and sample integrity are 3 important QC for each sample. Briefly describe how they can be done. **[1.5 marks]**

11. Why is validation and reproducibility so important in clinical pipelines? **[1 mark]**


## Section 2: Population Genomics

1. What command would you use to get the total number of variants from a VCF called file.vcf.gz? And then only for the region 10,000,000-30,000,000 from chromosome 3?

2. What information would the command bcftools query -l file.vcf.gz give you? (Look at the bcftools documentation for a complete list of options)

3. If a VCF file has GT:GQ:DP:GL in the FORMAT field for a particular variant, what kind of information about the samples is recorded?

4. Which of the following programs can convert a VCF file into data formats compatible with population genomics analyses?
	- a. PLINK
	- b. ADMIXTOOLS
	- c. CONVERTF
	- d. SMARTPCA

5. What was the extension of the three output files when you converted a VCF file using PLINK? Describe briefly what information is recorded in each file.

6. When using PLINK, what would the option --maf 0.10 do? (Look at the PLINK documentation for a complete list of options... and use the search box!)

7. The three EIGENSTRAT output files generated when converting VCF files with CONVERTF contain information about the genotypes, variants, and samples. How is the genotype information coded? Briefly describe each code key (4 in total).

8. Data missingness is characteristic of ancient DNA datasets and could lead to biases when building a PCA. If you use SMARTPCA on a dataset that includes both contemporary (no missing data) and ancient genotypes (missing data), what options can you use to avoid biases due to missing data?

9. When using population genomics data, PCA (select the correct answer):
	- a. is a formal test of shared ancestry between populations
	- b. should only be used as an exploratory tool to visualise genetic diversity and formulate hypotheses about population ancestry
	- c. a and b are correct
	- d. a and b are wrong

10. F and D statistics allow us to test hypotheses about populations ancestry:
	- a. True
	- b. False

11. You performed a F4 or a D statistic and the result is not significantly different from 0. What does it mean in terms of genetic admixture?


## Section 3: Structural Variation - **? marks**

Data for this part is located in the `/data/assignment3/` directory

## Section 3 Part 1: Copy number variations  (10 marks)


**Q3.1a** (2 marks)

![Copy Number](../images/Q1_fig1.tumour_CN.png)

The above figure shows the copy number estimates for various genomic regions inside a hypothetical human tumour sample.

List all the copy number variations that you can see from the figure.

**Q3.1b** (2 marks)

![VAF](../images/Q1_fig2.tumour_VAF.png)

Similarly, the above figure shows the variant allele frequency of variants inside another hypothetical tumour sample. 

List all the CNVs that you can see from this figure.


**Q3.1c** (2 marks)

Now interpret the two figure together as data from the same sample.

List the CNVs that you can see by integrating both data sets.


**Q3.1d** (2 marks)

The above figures are from an hypothetical "pure" tumour sample. But in practice, we often get tumour samples which are mixed with some normal/non-tumourous cells.

For example, for the same tumour sample, now there is a percentage of normal cells included in the data, and the CN estimates and VAF graphs now look like this:

![impure CN](../images/Q1_fig3.impure_CN.png)
![impure VAF](../images/Q1_fig4.impure_VAF.png)

Estimate tumour purity of the sample from these figures.

**Q3.1e** What do the grey bars in the figures represent? Why do they not contain any data points? - 1 mark
#=============== The below could be worth only 1, unclear from assignment

**Q3.1f.** If we suspect that p-arm of chromosome 8 has fused with the q-arm of chromosome 17. How will you be able to determine whether a fusion event has occured? Will short-read NGS data (WGS, WES or RNA-seq) be able to detect such a fusion event? - 2 marks


## Section 3 Part 2 - (10 marks)

**Q3.2a - data processing (6 points)**

In `~/data/assignment3/Q2` you should find these files:

```
sample_R1.fastq.gz
sample_R2.fastq.gz
mappable_region.fasta
mappable_region.fasta.fai
```

Download these files and write a BASH script to do the following:

(Short-read alignment, refer to Wk 5 Prac)

1. Use `bwa` to index the `mappabale_region.fasta` file. This is your reference sequence.

2. Map the paired-end FASTQ set (`sample_R1.fastq.gz` and `sample_R2.fastq.gz`) to the indexed reference above. Remember the output should be in BAM format, and needs to be sorted and indexed.

(Structural variation detection, refer to Wk 9 Prac)

3. Use `manta` to call structural variation on this set of data.


Submit:
* the BASH script
* the mapped and sorted BAM file
* the output file (`diploidSV.vcf.gz`) generated by Manta


**Q3.2b - interpreting the results (4 points)**

Either examine the read data on IGV or interpret manta output:

1) Locate and list any breakpoints you can see in the data,

2) Identify all SV events and their associated breakpoints. Show the steps and reasonings for your answers. Include diagrams if you think it helps. If you want to use hand-drawn diagram, just take and submit a photo of your drawing, but make sure it's clearly legible.



# **Section 4: Eukaryotic genome assembly**


## Introduction

In this assignment, you aim is to de novo assemble one of the SMALLEST eukaryotic genomes, Encephalitozoon intestinalis. E. intestinalis belongs to Microsporidia, and it's a parasite (microbial fungi), which causes microsporidiosis (an oppotunistic intestinal infection that causes diarrhea and wasting in immunocompromised individuals, such as HIV). If you want to understand more about E. intestinalis, please hava a read at [wikipedia](https://en.wikipedia.org/wiki/Encephalitozoon_intestinalis). 

Although the genome of E. intestinalis is very small (~2.5 Mb), it has 11 chromsomes. If you want to know more about the genome statistics of E. intestinalis, please have a look at [here](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_024399295.1/).

## Data

You are provided with 8 fastq files containing sequencing reads using different sequencing platforms. These fastq files can be found in `~/data/assignment2/raw_data` folder, and are including:

| File(s)                                    | Platform | Coverage | Description                                            |
|--------------------------------------------|----------|----------|--------------------------------------------------------|
| illumina_SR_30x_1.fq, illumina_SR_30x_2.fq | Illumina | ~30x     | Paried-end short reads from Illumina Miniseq           |
| nanopore_LR_15x.fq                         | Nanopore | ~15x     | Long reads from Nanopore MinION                        |
| nanopore_LR_15x_filt.fq                    | Nanopore | ~15x     | Long reads (reads >= 10 kb) from Nanopore MinION       |
| nanopore_LR_30x.fq                         | Nanopore | ~30x     | Long reads from Nanopore MinION                        |
| pacbio_LR_15x.fq                           | PacBio   | ~15x     | Long reads from PacBio_SMRT Sequel II                  |
| pacbio_LR_15x_filt.fq                      | PacBio   | ~15x     | Long reads (reads >= 10 kb) from PacBio_SMRT Sequel II |
| pacbio_LR_30x.fq                           | PacBio   | ~30x     | Long reads from PacBio_SMRT Sequel II                  |

The original dataset can be found [here](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=594722).

The illumina dataset will be used to do genome survey analysis (genome size estimation), and then you will be generating an assembly from each of six long reads datasets using Flye (v2.8.1) and then comparing the quality of these assemblies.

You are also provided with the E. intestinalis reference (taken from [here](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_024399295.1/) if you want to have a look). 

It is in the `~/data/assignment2/DB` directory and is called `GCA_024399295.1_ASM2439929v1_genomic.fna`. 

In addition to the sequencing files, you will be also given two scripts which will be used to do sequence/genome statistics and genome survey analysis. These two scripts can be found in folder `~/data/assignment2/bin`, and are including:

| Script       | Description                                 | Link                                               |
|----------------|---------------------------------------------|----------------------------------------------------|
| assembly-stats | A light tool to do basic genome statistics  | https://github.com/sanger-pathogens/assembly-stats |
| genomescope.R  | A R script to do genome survey analysis     | https://github.com/schatzlab/genomescope           |

### Part 1, QC

**1.1** Run `assembly-stats` on each LR (Long Reads) dataset and find out the following info for each dataset:

* total number of bases *- 1 marks*
* number of reads *- 1 marks*
* average read length*- 1 marks*
* largest read length*- 1 marks*

**1.2** Predict which datasets will produce the best and worst assemblies. 
Don't worry if your predictions don't match up with your results later. 
Just try to justify your predictions based on the information you've collected and your current knowledge. *- 2 marks*

### Part 2

**2.1** Use `jellyfish` and `genomescope.R` to perform genome survey analysis. What is the estimated genome size? Provide the generated figure showing the fitted model for k-mer distribution (Hint: plot.png in your genomescope output folder). *- 4 marks*

### Part 3

**3.1** Write a bash script to assemble all 6 LR datasets using Flye. [Hint: Assembly all 6 LR datasets will take ~50 mins in total, so be patient if you see the Flye is running for a while.] *- 6 marks*

**3.2** Run `assembly-stats` on each assembled genome and find out the following info for the assembled genomes:

* draft genome size *- 1 marks*
* number of contigs *- 1 marks*
* largest contig length*- 1 marks*
* N50 *- 1 marks*

### Part 4

**4.1** Compare these assemblies with each other using QUAST and provide the report (Hint: report.pdf in your QUAST output folder). *- 2 marks*

**4.2** Comment on the contig length distribution. Is this what you expected? *- 2 marks*

**4.3** Explain why contiguity isn't a good measure of assembly accuracy but is still relevant to the overall assessment of assembly quality. *- 2 marks*

### Part 5

**5.1** Run BUSCO on your 6 assemblies using an appropriate lineage and create a comparison image using the `generate_plot.py` script. Note that this script will only work when the `BUSCO` conda environment is activated. Provide this image [Hint: the "busco_figure.png" file in short summaries folder]. *- 4 marks*

**5.2** Comment on the BUSCO results. Which assembly appears to be the best and which is the worst? *- 2 marks*

### Part 6

**6.1** Considering your findings so far, justify which assembly you think is the best and which is the worst.
If your findings don't match your predictions from Part 1, try to explain why this might be. *- 4 marks*

## Check-list for your assignment2 submission

* A document including answers to all questions
* A image/figure showing the k-mer distribution model fitting on genome survey analysis in Part 2.1
* A bash script to assemble all 6 LR datasets using Flye in Part 3.1
* A report from the QUAST assessment in Part 4.1
* A comparison figure showing the BUSCO scores for all 6 draft assemblies in Part 5.1.
