
# Eukaryotic Genome Assembly Part 1
{:.no_toc}

#### By Chelsea Matthews
{:.no_toc}

* TOC
{:toc}

# **1. Introduction/Background**

## 1.1 Eukaryotic Genome Assembly

Today and in the next practical we will be looking at *de-novo* assembly of eukaryotic genomes which are generally larger and more complex than prokaryotes.
The human genome falls into this category as it is a large diploid (~3.1Gbp haplotype), comprised of ~50% repetitive elements.

While the focus of this practical is technically on the assembly of large complex genomes, it's not feasible for us to actually assemble a large complex genome (i.e. human genome) because our VM's are far too small.
Assembling a human genome from an appropriate amount of long-read data takes hundreds to thousands of CPU hours and hundreds of Gb of RAM (see Flye assembly benchmarking [here](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-flye-benchmarks)).   
Therefore we'll be assembling a small eukaryotic genome instead - fission yeast (*Schizosaccharomyces pombe*). Even though this genome is fairly small (~12Mbp), it will still take about 20 minutes for the assembly to run,  even with very low (~5x) long read coverage.

Here is some useful information about the fission yeast genome:
- Reference genome: ASM294v2
- Number of chromosomes: 3 nucleus chromosomes
- Genome size (reference): 12,591,251 bp
- Ploidy: Haploid

## 1.2 Practical Overview

We have four long-read datasets available. They are Nanopore reads with 5x and 10x coverage and PacBio reads with 5x and 10x coverage. We will assemble one of these datasets  in class (because it takes about 20 minutes to run) and the remaining three will be provided. We will then compare the quality of these assemblies using two different measures.

The steps in this analysis and their corresponding subdirectory are shown in the table below: 

| Subdirectory | Step                                |
| ------------ | ----------------------------------- |
| 0_data       | Get data                            |
| 1_qc         | Assess read quality and trim reads  |
| 2_assemble   | Assemble the genome from long reads |
| 3_quast      | Assess assembly contiguity          |
| 4_busco      | Assess completeness of gene space   |

## 1.3 Learning Outcomes
- Learn how to do basic quality assessment on long reads
- Learn how to do *de-novo* genome assembly using Flye with long reads
- Learn how to assess genome assembly quality using QUAST and BUSCO

# **2. Setup**

Activate your `bioinf` environment.

```bash
source activate bioinf
```

Create the directory structure for the prac and move into it.

```bash
mkdir -p euk_assembly_pt1/{0_data,1_qc,2_assemble,3_quast,4_busco}
cd euk_assembly_pt1
```

Now create symlinks to the data for todays practical in your `0_data` directory.

```bash
ln -s /shared/data/euk_assembly/part1/*.fq 0_data/.
# check contents of 0_data
ls -lh 0_data/
```

You should have symlinks to the four fastq files described in the table below. 

| File(s)         | Platform | Coverage | Description                         |
| --------------- | -------- | -------- | ----------------------------------- |
| nanopore_5x.fq  | Nanopore | ~5x      | Long reads from Nanopore PromethION |
| nanopore_10x.fq | Nanopore | ~10x     | Long reads from Nanopore PromethION |
| pacbio_5x.fq    | PacBio   | ~5x      | Long reads from PacBio_SMRT Sequel  |
| pacbio_10x.fq   | PacBio   | ~10x     | Long reads from PacBio_SMRT Sequel  |

We will also be using `chopper` today which isn't installed in the `bioinf` environment. 
Create a `scripts` directory and copy the `chopper-linux` script into it. Then change the script permissions (by giving all users executable rights) so that you can run it. 

```bash
mkdir scripts
# get the chopper script
cp /shared/data/euk_assembly/scripts/chopper-linux scripts/.
# make it executable
chmod +x scripts/chopper-linux
# check that it is executable
ls -lh scripts/chopper-linux
```


# **3. QC and Trimming**

## 3.1 Data quantity and quality
We'll first use `seqkit` to get some basic information about the reads in each of our `fastq` files. 

```bash
## First, make a file and add a header to it
echo "# Raw data" > 1_qc/summary.stats

## Run seqkit stats and append the results to the file we just made
seqkit stats 0_data/nanopore_5x.fq 0_data/nanopore_10x.fq 0_data/pacbio_5x.fq 0_data/pacbio_10x.fq -N 10,50,90 >> 1_qc/summary.stats

```

Use `less`  to view the file you created to answer the questions below. You will need to do some basic calculations. 

❓**Questions:**
- Which dataset has the largest average read length?
* Which dataset has the longest individual read?
* How would you calculate the read coverage for a single dataset (assuming a 12.5Mbp genome)?
* What information does the N90 column provide?

That was a very high level summary of our data but it would be good to be able to check the quality of our reads. 

We will only run these steps for the two 10x datasets to save time. In any other scenario, we would assess the quality of all of our data. 

We will be using the tool `NanoPlot` for quality assessment. Despite its name, it works for PacBio reads too. We will only be running it on the two 10x datasets to save time. However, in any other scenario, you would be assessing the quality of all of your data. 

Run `NanoPlot` on  your 10x long-read datasets as below.  The directories are provided for all four datasets but the code for the 5x samples is commented out.   

```bash
mkdir -p 1_qc/{nanopore_5x,nanopore_10x,pacbio_5x,pacbio_10x}

#NanoPlot -t 2 --outdir 1_qc/nanopore_5x --prefix nanopore_ --fastq 0_data/nanopore_5x.fq

NanoPlot -t 2 --outdir 1_qc/nanopore_10x --prefix nanopore_ --fastq 0_data/nanopore_10x.fq

#NanoPlot -t 2 --outdir 1_qc/pacbio_5x --prefix pacbio_ --fastq 0_data/pacbio_5x.fq

NanoPlot -t 2 --outdir 1_qc/pacbio_10x --prefix pacbio_ --fastq 0_data/pacbio_10x.fq

```

NanoPlot produces quite a lot of output and a report that is (sort of) similar to the report produced by FASTQC. Use the file browser to find and open both the `nanopore_NanoPlot-report.html`  and `pacbio_NanoPlot-report.html`  files in a web browser. 
Try to answer the questions below.

❓**Questions:**
- Compare the "Non weighted histogram of read lengths" plots (the third one down) between your Nanopore and PacBio datasets. What do you see? 
 - Looking at the very last plot in the report - "Read lengths vs Average read quality kde plot", in what range do most of the average read qualities for Nanopore reads fall? What about PacBio reads?
 - Do you think that removing all reads with average quality less than 10 would be appropriate for the Nanopore dataset? Discuss.

Keep in mind that this particular Nanopore dataset is quite old and reads generated more recently is generally of higher quality than this. PacBio read quality has also improved and is still more accurate than Nanopore but the difference between the two technologies is not as big as we see in these older datasets. 

## 3.2. Trimming

Run `chopper` as below to remove any reads with an average quality of less than 10. Nothing will be removed from the PacBio datasets but this will ensure all of your trimmed data is in the same location.

```bash
./scripts/chopper-linux -q 10 -i 0_data/nanopore_5x.fq > 1_qc/nanopore_5x.fq

./scripts/chopper-linux -q 10 -i 0_data/nanopore_10x.fq > 1_qc/nanopore_10x.fq

./scripts/chopper-linux -q 10 -i 0_data/pacbio_5x.fq > 1_qc/pacbio_5x.fq

./scripts/chopper-linux -q 10 -i 0_data/pacbio_10x.fq > 1_qc/pacbio_10x.fq
```

Now, run `seqkit` on your trimmed fastq files again to see how much sequencing data (and hence approximate genome coverage) you have remaining. 

```bash
echo -e "#\n# Trimmed data" >> 1_qc/summary.stats

seqkit stats 1_qc/nanopore_5x.fq 1_qc/nanopore_10x.fq 1_qc/pacbio_5x.fq 1_qc/pacbio_10x.fq -N 10,50,90 >> 1_qc/summary.stats 
```

❓**Question:**
- How have the avg_len and max_len for the Nanopore reads changed? 
- How much read coverage do you have now in your Nanopore 5x dataset (assuming a 12.5Mbp genome). 
- Do you think that you have enough data to generate a good quality assembly?
# **4. Assembly**

Now we're ready to assemble our genome. 
We will only be running one assembly during this practical as they take at least 20 minutes to run. 

There are many different assembly tools available which include:

* Canu 
* Flye
* Hifiasm
* NextDenovo
* Wtdbg2 
* Raven
* Falcon
* Shasta
* miniasm

These tools all do more or less the same job (ie. long read assembly) but with slight differences. 

## 4.1 Assemble Genome

We will be using Flye for our assembly. 
Check that it is working as below. The usage instructions should be printed to the screen.

```bash
flye
```

We will be assembling just one of our datasets in class and assemblies generated from the remaining three datasets will be provided. Code is provided for the nanopore_5x dataset. 

```bash
mkdir -p 2_assemble/nanopore_5x

flye --nano-raw 1_qc/nanopore_5x.fq --out-dir 2_assembly/nanopore_5x --threads 2
```

This will take about 20 minutes to run. 

### While we wait...

One of the difficulties with assembling larger genomes is working out what resources a tool will require and deciding whether it will be suited to assembling your genome of interest within the resources available to you. 
The authors of Flye did some benchmarks on computational resources required when assembling genomes from different species with different input data. You can check this from this [link](https://github.com/fenderglass/Flye#flye-benchmarks). Answer the following questions by looking at the benchmarking table.

❓**Questions:**
* How many CPU hours (approximately) are required to assemble a bacteria?
* How many CPU hours (approximately) are required to assemble a mammal with high raw-read coverage?
* How many CPU hours (approximately) are required to assemble a mammal from HiFi reads?
* Why is the number of CPU hours and memory requirement lower for a HiFi assembly compared with a non-HiFi read assembly?
* If you had access to a single node with 48 compute threads and sufficient memory, how long (in hours) would it take to assemble a "well-behaved large genome"?
* How long would this same assembly take on your local VM (assuming memory was not a limitation) with two compute threads?


## 4.2 When Flye finishes...

When flye finishes, it will  print some basic information about your assembly to the terminal. 

❓**Questions:**
- How long is your assembly in total? 
- How many contigs/fragments are there?
- If your assembly was perfect, how many fragments would you expect?

Flye will also generate 5 folders, which store output files from the 5 stages of Flye, and some important individual files. These include:

- **assembly.fasta**: This is the final assembly. Contains contigs and possibly scaffolds.
- params.json: the parameters that Flye used for this run
- assembly_graph.{gv or gfa}: Final repeat graph. The edges of repeat graph represent genomic sequence, and nodes define the junctions. You can get more info about the repeat graph from [here](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-repeat-graph) 
- assembly_info.txt: Extra information about contigs.
- flye.log: Log report showing all running info and summary of final assembly.


# **5. Assembly Quality - QUAST**

Now that we have an assembly, let's look at assessing assembly quality.

The tool QUAST (QUality ASsesment Tool) [documentation here](http://quast.sourceforge.net/docs/manual.html) can be run on one or more assemblies at once and produces some handy comparison statistics and graphics.
It can be run with or without a reference genome. 
If we provide a reference genome, each of the assemblies will be compared to this reference which allows QUAST to produce some additional statistics. 

Here we will run it without a reference.

Because you have only generated one assembly, I've provided four other assemblies so that you have something to compare with. Place all of these assemblies in a single directory `2_assembly/all` so that they're easy to find. Both of the provided Nanopore assemblies (5x and 10x) were generated from the untrimmed datasets.

```bash
mkdir -p 2_assembly/all

# copy your assembly to the new dir and rename
cp 2_assembly/nanopore_5x/assembly.fasta 2_assembly/all/nanopore_5x.fasta

# Copy provided assemblies to new dir
cp /shared/data/euk_assembly/part1/assemblies/* 2_assembly/all/.
```

Now we'll run QUAST on all 5 assemblies.

```bash
quast -o 3_quast/without_ref -t 2 --labels "nanopore_5x,nanopore_10x,pacbio_5x,pacbio_10x" 2_assembly/all/nanopore_5x.fasta 2_assembly/all/nanopore_10x.assembly.fasta 2_assembly/all/pacbio_5x.assembly.fasta 2_assembly/all/pacbio_10x.assembly.fasta

```

This won't take very long (< 1 min).

Open the `icarus.html` file located in the `3_quast/without_ref` directory in a web browser. 

Click on the QUAST Report link and have a look at the metrics produced by QUAST for each of your assemblies. 

* Based on the N50 and N90 scores, which is the best assembly?

You've previously learnt about N50 as a measure of assembly contiguity ([refresh your memory here](https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics#N50)) but N50 on its own only gives us a snapshot into assembly contiguity.
A better way to look at assembly contiguity is to inspect a cumulative contig length plot.
This is particularly helpful when we want to compare multiple different assemblies of the same genome or species. 

Scroll down and inspect the Cumulative Length plot. Look at the Nx plot and GC Content plots too while you're there. 

❓**Questions:**
* Given that the fission yeast genome should be about 13Mbp long, which assembly do you think looks the best based on the contig length distributions?
* Spend some time interpreting the cumulative length plot. Under what circumstances would a cumulative length plot be more informative than just comparing N50s?
* What would an assembly with only one contig look like on this plot?

Now let's have a look at the Icarus contig size viewer. 
Click on the "View in Icarus contig browser" link at the top of the page. 

This isn't really new information but it can be nice to visually see the lengths of our contigs. 
It's also a nice way to see the N50 and N90 values distributed on a visual representation of the assembly contigs. 

* Based on what you've seen in the QUAST report, which assembly do you think looks the best?

# **6. Assembly Quality - BUSCO**

Firstly, BUSCO is not installed in the `bioinf` environment. It is in an environment called `busco`. Therefore, to use it, we first have to activate the `busco` environment. This will make all of the software in the `bioinf` environment unavailable until you activate it again. 

```bash
source activate busco
```

BUSCO is a tool that helps us measure how well we have assembled the gene space of an assembly. 
It does this by searching our assembly for a list of genes that should be present in our assembly in a single copy. 
These genes are known as "Universal Single Copy Orthologs" and are genes that are present in at least 90% of the species/clade members and are only present in a single copy within at least 90% of the species/clade.
It is essentially a list of genes that we are almost certain should be present within our genome once.

Obviously this list of genes is different for different species/clades and so BUSCO has a list of datasets that we can choose from. 

Have a look at this list using the command below. 

```
busco --list-datasets
```

Before we run BUSCO, we need to decide which lineage dataset we should use for our assembled species. Get some ideas about the taxonomy of your assembled species is always a good starting point. For example, for fission yeast, we can get its taxonomy from [NCBI taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4896). 

Compare the taxonomy with the above listed available BUSCO lineage datasets and you'll find that there are three datasets we should be able to use for fission yeast, which are `eukaryota_odb12.2`, `fungi_odb12.2` and `ascomycota_odb12.2` following the taxonomy tree.
Datasets closer to the branch end of taxonomy will containsingle copy orthologs that are more specific to the target species.

❓Why does the `ascomycota` database have more single-copy orthologs than the `eukaryota` database? 

Generally speaking, you should choose the lineage dataset closer to the branch end of taxonomy. But in this prac, choosing `ascomycota` will take much longer time to run than choosing `eukaryota` because there are a so many more orthologs. 
Therefore, we will use `eukaryota`.

**NOTE:** The databases listed by BUSCO are mostly `odb12.2` because that is the most recent version available online. HOwever, our installation includes only `odb10`. Therefore,  make sure that you specify `eukaryota_odb10` when you run BUSCO. 

Let's run BUSCO on your Nanopore 5x assembly. 

```bash
busco -m genome -i 2_assembly/nanopore_5x/assembly.fasta -o 4_busco/nanopore_5x -f -c 2 -l eukaryota_odb10
```

This will take ~10 mins to finish. 
When the job is finished, you should get output printing to the terminal that includes something similar to the following. Let's look at the example here while we wait.

```
        --------------------------------------------------
        |Results from dataset eukaryota_odb10             |
        --------------------------------------------------
        |C:94.9%[S:92.9%,D:2.0%],F:2.0%,M:3.1%,n:255      |
        |242    Complete BUSCOs (C)                       |
        |237    Complete and single-copy BUSCOs (S)       |
        |5      Complete and duplicated BUSCOs (D)        |
        |5      Fragmented BUSCOs (F)                     |
        |8      Missing BUSCOs (M)                        |
        |255    Total BUSCO groups searched               |
        --------------------------------------------------
```

These results are also stored in the `short_summary*.txt` file in `4_busco/nanopore_5x` directory and are what we need for the next step. 

Now, we could interpret these results on their own but it's a lot more interesting if we have more than one assembly to compare. 

I have already run BUSCO on each of the provided assemblies. Make a new directory, copy over your assembly results, and then the BUSCO results for the other three assemblies. We put them all in the same directory because the BUSCO `generate_plot.py` script (this is usually installed along with BUSCO) requires it. 

```bash
mkdir -p 4_busco/all 

# copy your busco results to the new dir
cp 4_busco/nanopore_5x/short_summary*.txt 4_busco/all/.

# copy provided busco results to new dir
cp /shared/data/euk_assembly/part1/busco_results/* 4_busco/all/.
```

Do a quick check in the `4_busco/all` folder to check that there are four files named like below:

`short_summary.specific.eukaryota_odb10.xxx.txt`

Now let's visualise the results. 

```bash
generate_plot.py -wd 4_busco/all
```

If all goes well, this should produce a `busco_figure.png` in the `4_busco/all` directory. 

Navigate to it and open the png. 

* Which assembly is the best according to the BUSCO metrics? Why?
* *S. pombe* has a haploid genome but there are duplicated genes present in the assembly. What do you think this means? 
* Based on the BUSCO and QUAST results, which assembly do you think looks the best?
- Do you think that the best of these assemblies would be considered a reference-quality assembly? Discuss.


