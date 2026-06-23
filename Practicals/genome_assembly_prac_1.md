# [Genomics](http://university-of-adelaide-bx-masters.github.io/Genomics/)
{:.no_toc}

# Genome Assembly Practical: Part 1
{:.no_toc}

*By Zhipeng Qu and Chelsea Matthews* 

* TOC
{:toc}

# **Introduction**

## 1.1 *De-novo* assembly of eukaryotic genomes

Today and in the next practical we will be looking at *de-novo* assembly of eukaryotic genomes which are generally larger and more complex than prokaryotes.  
The human genome falls into this category as it is a large diploid (~3.1Gbp haplotype), comprised of ~50% repetitive elements. 

## 1.2 Practical Overview

While the focus of this practical and the next is on the assembly of large complex genomes, it's not feasible for us to actually assemble a large complex genome (i.e. human genome) because frankly, our VM's are far too small.
Assembling a human genome from an appropriate amount of long-read data takes hundreds to thousands of CPU hours and hundreds of Gb of RAM (see Flye assembly benchmarking [here](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-flye-benchmarks)).   
Therefore we'll be assembling a small eukaryotic genome, fission yeast (*Schizosaccharomyces pombe*).
Even though the fission yeast genome is fairly small (~15 Mb with 3 chromosomes), it's still more complicated than the prokaryotic organisms you've looked at so far and it will still take ~20 minutes for the assembly to run even with very low (~5x) long read coverage.

What we will actually do: First estimate the size of a haploid yeast genome using illumina reads. 
Then Denovo assemble HIFi and ONT and check quality of resulting assemblies.
This will highlight the challenges caused by repetitive regions + genome size in comparison with long reads. 


## 1.3 Learning Outcomes 

- Practice bash commands learned previously
- Practice QC for Illumina data
- Learn how to do *de-novo* genome assembly using Flye with long reads
- Learn how to assess your genome assembly quality using QUAST and BUSCO
- Understand the importance of haplotype phasing and the concept of trio-sequencing for phasing haplotypes

# **2. Setup**

## 2.1 Software

As in previous weeks, you will be using RStudio to interact with your VM.
Let's activate the 'bioinf' environment so that we can access the software we'll need for the practical.

```bash
source activate bioinf
```

The table below lists all of the tools we will be using in this Practical and the next. 

| Tool/Package    | Version      | URL                                                        |
|----------------|--------------|------------------------------------------------------------|
| fastQC         | v0.11.9      | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| assembly-stats | v1.0.1       | https://github.com/sanger-pathogens/assembly-stats         |
| jellyfish      | v2.2.10      | https://github.com/gmarcais/Jellyfish                      |
| genomescope    | v1           | https://github.com/schatzlab/genomescope                   |
| flye           | v2.8.1-b1676 | https://github.com/fenderglass/Flye                        |
| QUAST          | V5.2.0       | https://github.com/ablab/quast                             |
| BUSCO          | v5.4.4       | https://busco.ezlab.org/                                   |

Two of these tools, assembly-stats and genomescope, are not installed in the `bioinf` environment and will instead be run directly from scripts. Setting this up is covered in the next section. 


The following table shows the estimated run time on our VMs for the different processes we'll be running. 


| Step              | Tool/Package        | Estimated run time  |
|-------------------|---------------------|---------------------|
| QC                | fastqc              | < 1 min             |
| QC                | assembly-stat       | < 1 min             |
| genome survey     | jellyfish           | < 5 mins            |
| genome survey     | genomescope         | < 1 min             |
| genome assembly   | flye + nanopore_5x  | ~20 mins            |
| genome assembly   | flye + pacbio_5x    | ~20 mins (optional) |
| genome assembly   | flye + nanopore_10x | ~30 mins (optional) |
| genome assembly   | flye + pacbio_10x   | ~30 mins (optional) |
| genome assessment | QUAST               | < 5 min             |
| genome assessment | BUSCO               | ~10-20 mins (each)  |



## 2.2 Create directory structure

The following is the directory structure we'll be working within today and in the next practical. 


Create these directories using the code below. 

```bash
cd ~/
mkdir prac_genome_assembly
cd prac_genome_assembly
mkdir 01_bin 02_DB 03_raw_data 04_results 05_scripts
cd 04_results
mkdir 01_QC 02_genome_survey 03_genome_assembly 04_genome_assessment

```

Check your folder structure:

```bash
cd ~/
tree ./prac_genome_assembly
```

Because `assembly-stats` and `genomescope` are not installed in the `bioinf` environment, we will run them directly from scripts. 
The scripts are located in `~/data/prac_genome_assembly/01_bin/`, and we will put them in the `01_bin` folder of our project.

```bash
cd ~/prac_genome_assembly/01_bin
cp ~/data/prac_genome_assembly/01_bin/* ./
```


## 2.3 Get data

Create symlinks to the raw data needed for this practical as below. 
The original dataset from which this data was subset can be found [here](https://www.ncbi.nlm.nih.gov/sra?term=SRP352919).


```bash
cd ~/prac_genome_assembly/02_DB
cp ~/data/prac_genome_assembly/02_DB/* ./
cd ~/prac_genome_assembly/03_raw_data
cp ~/data/prac_genome_assembly/03_raw_data/*.fq ./
```

We will be using raw sequencing data from different sequencing platforms in this Prac. 
These fastq files can be found in `~/data/prac_genome_assembly/01_raw_data` and are described in the table below:

| File(s)                                    | Platform | Coverage | Description                                                |
|--------------------------------------------|----------|----------|------------------------------------------------------------|
| illumina_SR_20x_1.fq, illumina_SR_20x_2.fq | Illumina | ~20x     | Paried-end (PE150) short reads from Illumina MGISEQ-2000RS |
| nanopore_LR_5x.fq                          | Nanopore | ~5x      | Long reads from Nanopore PromethION                        |
| nanopore_LR_10x.fq                         | Nanopore | ~10x     | Long reads from Nanopore PromethION                        |
| pacbio_LR_5x.fq                            | PacBio   | ~5x      | Long reads from PacBio_SMRT Sequel                         |
| pacbio_LR_10x.fq                           | PacBio   | ~10x     | Long reads from PacBio_SMRT Sequel                         |

We will use the paired-end illumina short reads to do genome survey analysis (estimate the genome size), and use the other four Long Reads (LR) files to do genome assembly separately.

Due to limitation of computing resources in our VMs, we will use subsets of raw sequencing reads. 
The original sequencing dataset is very big, you can access it from this [link](https://www.ncbi.nlm.nih.gov/sra?term=SRP352919). Here is some useful information about the fission yeast genome:

- Reference genome: ASM294v2
- Number of chromosomes: 3 nucleus chromosomes
- Genome size (reference): 12,591,251 bp
- Ploidy: Haploid



# **3. Genome Survey Analysis**

When we start a whole genome sequencing project for a new species, we normally need to collect some genomics information before we do the de novo genome assembly, such as we need to know how big the genome is. In the lecture, we had learned that we can use lab-based flow cytometry to estimate the genome size, and we can also use short reads (illumina reads) to do this computationally. In this part, we will learn how to estimate the genome size using short reads.


## 3.1 QC of Illumina reads

The first step in any bioinformatics analysis is always quality control. 
Let's check the quality of our short reads (illumina reads) using `fastQC`.

```bash
cd ~/prac_genome_assembly/04_results/01_QC
fastqc ~/prac_genome_assembly/03_raw_data/illumina_SR_20x_1.fq ~/prac_genome_assembly/03_raw_data/illumina_SR_20x_2.fq -o ./ -t 2
```

* *How many sequences are there in the dataset?*
* *How long are the Illumina reads?*
* *How good are the illumina reads?* 


## 3.2 Get k-mer distribution

The first step of genome size estimation is to get the k-mer distribution using the available short reads. We can use `jellyfish` to do this.

```bash
cd ~/prac_genome_assembly/04_results/02_genome_survey

jellyfish count -C -m 21 -s 4G -o illumina_SR_20x.21mer_out \
    ~/prac_genome_assembly/03_raw_data/illumina_SR_20x_1.fq \
    ~/prac_genome_assembly/03_raw_data/illumina_SR_20x_2.fq

jellyfish histo -o illumina_SR_20x.21mer_out.histo illumina_SR_20x.21mer_out

```

`jellyfish` will break short reads into fixed length short sequences, which we call them [k-mers](https://en.wikipedia.org/wiki/K-mer) (we use 21-mer in this project). The size of k-mers should be large enough allowing the k-mer to map uniquely to the genome (a concept used in designing primer/oligo length for PCR). However, too large k-mers leads to overuse of computational resources. `21` is normally a good start. In the `jellyfish count` command, `-C` means we count k-mers at both strands, `-s 4G` is used to control the memory usage, and `-m 21` means we will count 21-mers. After we count k-mers, we use `jellyfish histo` to get the frequency of k-mers with certain copy numbers.

## 3.3 Estimate genome length 

After we get the `histo` file from the `jellyfish` run, we can use that to do genome survey with `genomescope`. `genomescope` is a R script, we can run it with following command:

```bash
cd ~/prac_genome_assembly/04_results/02_genome_survey

Rscript ~/prac_genome_assembly/01_bin/genomescope.R illumina_SR_20x.21mer_out.histo 21 150 illumina_SR_20x.21mer
```

In the command, `21` means we are using 21-mer, `150` is the short reads length, and `illumina_SR_20x.21mer` will be the output folder. We can check the k-mer distribution by checking the file `plot.png`, which is normally located in the output folder `illumina_SR_20x.21mer`


# **4. Long Read summary statistics**

Let's use the tool `assembly-stats` to get more information about our long reads. 

```bash
cd ~/prac_genome_assembly/04_results/01_QC
~/prac_genome_assembly/01_bin/assembly-stats ~/prac_genome_assembly/03_raw_data/nanopore_LR_5x.fq
~/prac_genome_assembly/01_bin/assembly-stats ~/prac_genome_assembly/03_raw_data/pacbio_LR_5x.fq
~/prac_genome_assembly/01_bin/assembly-stats ~/prac_genome_assembly/03_raw_data/nanopore_LR_10x.fq
~/prac_genome_assembly/01_bin/assembly-stats ~/prac_genome_assembly/03_raw_data/pacbio_LR_10x.fq
```

From the output text, answer the following:

* *Which dataset has the largest average read length?*
* *Which dataset has the longest individual read?*
* *Assuming that our genome is approximately 15 Mb long, what coverage do these reads give us? Remember that Coverage = (total number of bases in reads)/genome size*
* *Do you think this is sufficient to generate a good assembly? Why or why not?*


# **5. *De-novo* genome assembly** 

## 5.1 assembly using Flye

Okey! Now we can do the de novo genome assembly. There are tons of assemblers to do de novo genome assembly. Some popular ones are `canu`, `Flye`, `smartdenovo`, `wtdbg`, etc. In this prac, we will be using `Flye`.

Do a test run to see if it works using the following (make sure you are under the `bioinf` conda environment):

```
flye
```

The Flye usage instructions should be printed to the screen. 

Now let's using the flye to do assembly on the `nanopore_LR_5x` dataset. 

```bash
cd ~/prac_genome_assembly/04_results/03_genome_assembly

flye --nano-raw ~/prac_genome_assembly/03_raw_data/nanopore_LR_5x.fq --out-dir nanopore_LR_5x --threads 2
```

It should take ~20 minutes for the assembly to complete.   

## 5.2. While we wait ...  Resources for genome assembly

Larger, more complex genomes contain high percentages of repeats. Assembling these repeats accurately with short reads is not effective. Long reads that are able to span these repeat regions result in much better assemblies. There are many long read assembly tools available including the following: 

* Canu 
* Flye
* NextDenovo
* Wtdbg2 
* Raven
* Falcon
* Shasta
* miniasm

These tools all do more or less the same job (ie. long read assembly) but with slight differences. 

When we are assembling a small genome, memory (RAM) and CPU (compute threads/cores) requirements are relatively low and so we can choose our assembler based almost entirely on the quality of the resulting assembly. 
While we can still do this to some degree for larger genomes, resource allocations and limitations in High Performance Computing (HPC) environments mean that some tools may not be suitable for our particular dataset (we may not be able to run the assembly to completion due to time or resource limitations) and we need to more carefully consider the parameters we use for the assembly so that we don't need to re-run it unnecessarily. 
Even without considering the cost of wasted resources (and this can be measured in dollars where we are using a service like Amazon Web Service (AWS)), re-running these assemblies can take a long time. 

One of the difficulties with assembling larger genomes is working out what resources a tool will require and deciding whether it will be suited to assembling your genome of interest within the resources available to you. 
The authors of Flye did some benchmarks on computational resources required when assembling genomes from different species with different input data. You can check this from this [link](https://github.com/fenderglass/Flye#flye-benchmarks). Check following questions after you have a look at the table.

* *How many CPU hours (approximately) are required to assemble a bacteria?*
* *How many CPU hours (approximately) are required to assemble a mammal with high raw-read coverage?*
* *How many CPU hours (approximately) are required to assemble a mammal from HiFi reads?*
* *Why is the number of CPU hours and memory requirement lower for a HiFi assembly compared with a non-HiFi read assembly?*
* *If you had access to a single node with 48 compute threads and sufficient memory, how long (in hours) would it take to assemble a "well-behaved large genome"?* 
* *How long would this same assembly take on your local VM (assuming memory was not a limitation) with two compute threads?*

You can see now why we aren't assembling a larger genome (i.e. human genome) using our VM. 

## 5.3 Looking at the assembly report

Flye will generate 5 folders, which store output files from 5 stages of Flye, and some important individual files. These includes:

- params.json: the parameters that Flye used for this run
- assembly_graph.{gv or gfa}: Final repeat graph. The edges of repeat graph represent genomic sequence, and nodes define the junctions. You can get more info about the repeat graph from [here](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-repeat-graph) 
- **assembly.fasta**: This is the final assembly. Contains contigs and possibly scaffolds.
- assembly_info.txt: Extra information about contigs.
- flye.log: Log report showing all running info and summary of final assembly.

You can check the log report using following commands:

```
cd ~/prac_genome_assembly/04_results/03_genome_assembly/nanopore_LR_5x
ls 

less flye.log
```
Type `G` to go to the end of the file, and you will see a brief summary about the assembly.

* *How many contigs does the final assembly have?*
* *How long is the final assembly and how does this compare with the estimated size of the genome?*

When you are finished, exit out of the report with `q`. 

### 5.4 Assembly using other datasets

You can do the genome assembly for the other three LR datasets:

```bash
cd ~/prac_genome_assembly/04_results/03_genome_assembly

flye --pacbio-raw ~/prac_genome_assembly/03_raw_data/pacbio_LR_5x.fq --out-dir pacbio_LR_5x --threads 2
flye --nano-raw ~/prac_genome_assembly/03_raw_data/nanopore_LR_10x.fq --out-dir nanopore_LR_10x --threads 2
flye --pacbio-raw ~/prac_genome_assembly/03_raw_data/pacbio_LR_10x.fq --out-dir pacbio_LR_10x --threads 2

```

Each of these will take 20-30 mins. I encourage you to run these jobs before the next session, and go through the reports (flye.log) to get ideas about the different assemblies, and we will look into them in more details in the next session.

### 5.5 Additional information - diploid assembly

The fission yeast is a small haploid genome which means that it is fairly straightforward nowadays to assemble. 
What is more challenging are larger diploid and polyploid genomes.

* *Can you think of a few reasons why diploid or polyploid genomes are more difficult to assemble than haploid genomes? Relate these reasons back to the OLC algorithm if you can.*
* *If we were assembling a diploid genome with a genome size of 2Gbp, under what circumstances might our resulting assembly be greater than 2Gbp in length?*

Canu is capable of assembling both haplotypes of a diploid genome - termed phasing. 
Have a look at the paper below (mainly the first figure) to understand how this works.

[Haplotype-Resolved Cattle Genomes Provide Insights Into Structural Variation and Adaptation](https://www.biorxiv.org/content/10.1101/720797v3.full)

