
# Eukaryotic Genome Assembly Part 2
{:.no_toc}

#### By Chelsea Matthews
{:.no_toc}

* TOC
{:toc}

# **1. Introduction/Background**

## 1.1 Eukaryotic Genome Assembly
In this practical we will be assembling a small portion of the New World screwworm - *Cochliomyia hominivorax* . The New World screwworm is a parasitic fly that lays eggs in open wounds or mucous membranes and when they hatch, the larvae burrow into healthy living tissue causing severe pain, tissue destruction, and may be fatal within a week if not treated. 

The New World screwworm is a **diploid** with a haploid genome size of ~534Mbp. We will use trio sequencing data to assemble and separate the maternal and paternal haplotypes in their offspring for a ~12Mbp region of chromosome 2. 

## Learning Outcomes

# **2. Setup**

Activate your `bioinf` environment.

```bash
source activate bioinf
```

Create the directory we'll be working in today and move into it.
```bash
cd ~
mkdir euk_assembly_pt2

cd euk_assembly_pt2
```

Create the directory structure shown below. 

```
euk_assembly_pt2/
├── 0_data
├── 1_trim
├── 2_assembly
├── 3_quast
└── 4_bandage

```

Copy all of the trio sequencing data located in `~/data/euk_assembly/part2/`into your `0_data` directory.

❓What files have you got?

We'll also be using a tool called `yak` that isn't installed in the `bioinf` environment. 
To use it, copy the `yak` directory into your local directory. 

```bash
cp -r /shared/data/yak .
```
# 3. Quality Control

## 3.1 PacBio

Run NanoPlot on your long reads as below and open the `NanoPlot-report.html` file in a web browser.

```bash
NanoPlot -t 2 --outdir 1_trim/hifi_30x --fastq 0_data/hifi_30x.fq 
```

❓**Questions:**
- What is the mean read length of your PacBio HiFi reads?
- How much coverage do you have of each haplotype, assuming that the genome (at least the portion we're working with) is ~12Mbp long.
- What does the "Weighted histogram of read lengths" plot tell us?
- Does there seem to be any relationship between read length and mean basecall quality?

We aren't going to trim our PacBio reads today because this data is very high quality but if we did want to, we could use `chopper` like in the previous practical.

## 3.2 Illumina

Let's run `fastqc` on our Illumina reads to see what we're working with. 

```bash
mkdir -p 0_data/fastqc/

fastqc 0_data/*.fq.gz -o 0_data/fastqc/ -t 2
```

Focusing on just the fastqc report for the "mother" reads, answer the following questions:

❓**Questions:**
- What is the length of these reads? What does this tell you?
- Do these reads have the characteristic 3' quality drop-off commonly seen in Illumina sequencing?
- If we trim our reads using a sliding window from the 3' end, what quality cutoff do you think might be appropriate?
- The reads analysed here are only a small portion of the reads generated in this sequencing run. Does this impact your interpretation of the "Per tile sequence quality" plot?
- Are these reads contaminated by adapters? 

Now let's trim. We'll use `fastp`. 
Adapter trimming is enabled by default in `fastp`.
Make sure to replace the `??`  with the quality cutoff you chose above for the `--cut_tail_mean_quality` option. 

```bash
fastp --in1 0_data/father_R1.fq.gz --in2 0_data/father_R2.fq.gz --out1 1_trim/father_R1.fq.gz --out2 1_trim/father_R2.fq.gz --cut_tail --cut_tail_mean_quality ?? --html 1_trim/father.html --json 1_trim/father.json

fastp --in1 0_data/mother_R1.fq.gz --in2 0_data/mother_R2.fq.gz --out1 1_trim/mother_R1.fq.gz --out2 1_trim/mother_R2.fq.gz --cut_tail --cut_tail_mean_quality ?? --html 1_trim/mother.html --json 1_trim/mother.json
```

Now, run `fastqc` again and check the report to make sure trimming worked as you expected. 

```bash
mkdir -p 1_trim/fastqc/
## fastqc trimmed reads
fastqc 1_trim/*.fq.gz -o 1_trim/fastqc/ -t 2
```

# **4. Assembly with Hifiasm**

Today we'll be using the genome assembly tool [Hifiasm](https://github.com/chhylp123/hifiasm). Hifiasm was originally designed specifically for assembling PacBio HiFi reads but now also supports Oxford Nanopore reads. It is one of the best haplotype-resolved assemblers for the trio-binning approach (using parental short reads) which is what we will be doing. 

Before running the actual assembly step with Hifiasm and our long reads, we need to run Yak to count k-mers in our short Illumina reads. 

```bash
# Run yak
./yak/yak count -b32 -t2 -o 2_assemble/father.yak <(zcat 1_trim/father_R1.fq.gz 1_trim/father_R2.fq.gz) <(zcat 1_trim/father_R1.fq.gz 1_trim/father_R2.fq.gz)

./yak/yak count -b32 -t2 -o 2_assemble/mother.yak <(zcat 1_trim/mother_R1.fq.gz 1_trim/mother_R2.fq.gz) <(zcat 1_trim/mother_R1.fq.gz 1_trim/mother_R2.fq.gz)
```

❓**Questions while yak runs:**
- What is a haplotype resolved assembly?
- Why are we counting k-mers in our Illumina reads?
- Why do we care about generating a haplotype resolved assembly?

Once `yak` has finished, we're ready to run `hifiasm`. 

```bash
mkdir -p 2_assemble/trio

hifiasm -f0 -o 2_assemble/trio -t2 -1 2_assemble/father.yak -2 2_assemble/mother.yak 0_data/hifi_30x.fq
```

This will take about 10 minutes to run. 

In the meantime, go to the [Hifiasm github](https://github.com/chhylp123/hifiasm#getting-started) and answer the following questions:

❓**Questions:**
- What does the `-f0` parameter in the command above do? (Ctrl + F in the github repo for `-f0`)
- Why do you think we're using this particular setting? 
- If we had Oxford Nanopore reads, what flag would we use in our `hifiasm` command?
- The Introduction in the Hifiasm github states "Hifiasm produces arguably the best single-sample telomere-to-telomere assemblies combing HiFi, ultralong and Hi-C reads". What properties do these data types have and how might they be used to generate such a high quality assembly?
- We are assembling a 12Mbp section of the New World Screwworm genome. Would it be appropriate to run BUSCO on our assembly?

Hopefully your assembly has completed! 

Let's check out what `hifiasm` has produced. 

```bash
ls -l1 2_assemble/trio/
```

Unlike `flye`, `hifiasm` doesn't generate a simple FASTA file with our assembly. Instead, it produces files in the `.gfa` format, standing for Graphical Fragment Assembly. 

Genome assemblers use graph structures (ie. a set of nodes connected by edges) to assemble genomes and these graph structures are often communicated using the .gfa format. We'll use `bandage` near the end of this practical to visualise these graph structures but for now, we can use the code below to extract just the genomic sequences of contigs from a `.gfa` file.  Can you tell what the command is doing? 

```bash
awk '/^S/{print ">"$2;print $3}' afile.gfa > output.fa
```

The files we're interested in are:

- `dip.hap1.p_ctg.gfa <- fully phased paternal haplotype contigs
- `dip.hap2.p_ctg.gfa` <- fully phased maternal haplotype contigs
- `prefix.p_ctg.gfa` <- assembly graph of primary contigs


Before we move on, let's generate another assembly with Hifiasm, this time with just our Hifi reads and without using the trio-binning approach. 

```bash
mkdir -p 2_assemble/just_hifi

hifiasm -f0 -o 2_assemble/just_hifi -t2 0_data/hifi_30x.fq
```
# **5. QUAST**

Now let's run QUAST on our three 


# **6. Bandage**

```bash
Bandage image 30x.dip.p_utg.noseq.gfa primary_unitigs.png
```


