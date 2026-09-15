
# Eukaryotic Genome Assembly Part 2
{:.no_toc}

#### By Chelsea Matthews
{:.no_toc}

* TOC
{:toc}

# **1. Introduction/Background**

## Eukaryotic Genome Assembly
In this practical we will be assembling a small portion of the New World screwworm - *Cochliomyia hominivorax* . The New World screwworm is a parasitic fly that lays eggs in open wounds or mucous membranes and when they hatch, the larvae burrow into healthy living tissue causing severe pain, tissue destruction, and may be fatal within a week if not treated. 

The New World screwworm is a **diploid** with a haploid genome size of ~534Mbp. We will use trio sequencing data to assemble and separate the maternal and paternal haplotypes in their offspring. However, the full genome will take far too long on our VMs so our data has been subset to a ~7Mbp region of chromosome 2. 

If you're interested, the paper associated with the original data is [here](https://academic.oup.com/g3journal/article/16/5/jkag053/8502022). 


## Learning Outcomes
- Practice data quality control 
- Learn how to run a simple trio-binning assembly pipeline
- Understand the differences between different types of phased assemblies (ie, collapsed, partially phased, haplotype-resolved)

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

```bash
 mkdir -p {0_data,1_trim,2_assembly,3_quast}
```

```
euk_assembly_pt2/
├── 0_data
├── 1_trim
├── 2_assembly
└── 3_quast

```

Copy all of the trio sequencing data located in `~/data/euk_assembly/part2/`into your `0_data` directory.

```bash
ls -lh 0_data
```

You should have 6 files. 
- Paired end Illumina reads from the mother and father
- 20x PacBio HiFi reads from their offspring
- 1 reference genome

We'll also be using a tool called `yak` (for k-mer analysis) that isn't installed in the `bioinf` environment. 
To use it, copy the `yak` directory into your local directory. 

```bash
cp -r /shared/data/yak .
```


# **3. Quality Control**

## PacBio

Run NanoPlot on your long reads as below and open the `NanoPlot-report.html` file in a web browser.

```bash
NanoPlot -t 2 --outdir 1_trim/hifi_30x --fastq 0_data/hifi_20x.fq.gz 
```

❓**Questions:**
- What is the mean read length of your PacBio HiFi reads?
- How much coverage do you have of each haplotype, assuming that the genome (at least the portion we're working with) is ~7Mbp long.
- What does the "Weighted histogram of read lengths" plot tell us?
- Does there seem to be any relationship between read length and mean basecall quality?

We aren't going to trim our PacBio reads today because this data is very high quality but if we did want to, we could use `chopper` like in the previous practical.

## Illumina

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

## Trio Assembly
Before running the actual assembly step with Hifiasm and our long reads, we need to run Yak to count k-mers in our short Illumina reads. 

```bash
# Run yak
./yak/yak count -b32 -t2 -o 2_assembly/father.yak <(zcat 1_trim/father_R1.fq.gz 1_trim/father_R2.fq.gz) <(zcat 1_trim/father_R1.fq.gz 1_trim/father_R2.fq.gz)

./yak/yak count -b32 -t2 -o 2_assembly/mother.yak <(zcat 1_trim/mother_R1.fq.gz 1_trim/mother_R2.fq.gz) <(zcat 1_trim/mother_R1.fq.gz 1_trim/mother_R2.fq.gz)
```

❓**Questions while yak runs:**
- What is a haplotype resolved assembly?
- Why are we counting k-mers in our Illumina reads?

Once `yak` has finished, we're ready to run `hifiasm`. 

```bash
hifiasm -f0 -o 2_assembly/trio -t2 -1 2_assembly/father.yak -2 2_assembly/mother.yak 0_data/hifi_20x.fq.gz

```

This will take only couple of minutes to run. 

In the meantime, discuss the following questions:

❓**Questions:**
- The Introduction in the Hifiasm github states "Hifiasm produces arguably the best single-sample telomere-to-telomere assemblies combing HiFi, ultralong and Hi-C reads". What properties do these data types have and how might they be used to generate such a high quality assembly?
- We are assembling a 7Mbp section of the New World Screwworm genome. Would it be appropriate to run BUSCO on our assembly?

## Assembly without trio data

Let's also assemble our genome from just our hifi reads, excluding the parental Illumina reads. 

```bash
hifiasm -f0 -o 2_assembly/hifionly -t2 0_data/hifi_20x.fq.gz
```

- What is the differences between these two approaches? How will this impact the resulting assemblies? ie. What type of assembly will each of these produce? See [Heng Li's blog](https://lh3.github.io/2021/04/17/concepts-in-phased-assemblies) for a refresher on the different types of assemblies. 
# **5. Hifiasm Output**

Let's check out what files `hifiasm` has produced. There's a lot. We want to focus mainly on the `.gfa` files and let's look at the `trio` assembly files first.

```bash
ls -lh 2_assembly
# print just the trio gfa files
ls -1 2_assembly/trio*.gfa
```

You should see something like this:
```txt
trio.dip.hap1.p_ctg.gfa
trio.dip.hap1.p_ctg.noseq.gfa
trio.dip.hap2.p_ctg.gfa
trio.dip.hap2.p_ctg.noseq.gfa
trio.dip.p_utg.gfa
trio.dip.p_utg.noseq.gfa
trio.dip.r_utg.gfa
trio.dip.r_utg.noseq.gfa
```

Unlike `flye`, `hifiasm` doesn't generate a simple FASTA file of our assembly. Instead, it produces files in the `.gfa` format, standing for Graphical Fragment Assembly. If you look carefully, the `flye` assembler also generates an `assembly.gfa` file, we just didn't look at it. 

Genome assemblers build graph structures (ie. a set of nodes connected by edges) from the provided reads to reconstruct the genome and these graph structures are often communicated using the `.gfa` format.

We are interested mainly in the `hap1` and `hap2`  graphs. And you can see that each one has a corresponding `*noseq.gfa` which is the same graph structure but with the genomic sequence removed to enable fast visualisation. 

The files we're most interested in are:
- `trio.dip.hap1.p_ctg.gfa` <- fully phased paternal haplotype contig graph
- `trio.dip.hap2.p_ctg.gfa`  <- fully phased maternal haplotype contig graph

To get a `fasta` file for each of the two haplotypes, use the code below. 

```bash
cd 2_assembly

awk '/^S/{print ">"$2;print $3}' trio.dip.hap1.p_ctg.gfa > trio.hap1.fa

awk '/^S/{print ">"$2;print $3}' trio.dip.hap2.p_ctg.gfa > trio.hap2.fa
```

Now, let's do the same for our `hifionly` assembly. 

```bash
# take a look first
ls -lh hifionly*.gfa

# get hap1
awk '/^S/{print ">"$2;print $3}' hifionly.bp.hap1.p_ctg.gfa > hifionly.hap1.fa

# get hap2
awk '/^S/{print ">"$2;print $3}' hifionly.bp.hap2.p_ctg.gfa > hifionly.hap2.fa

cd ..
```

# **6. QUAST**

Let's run QUAST on these files to get some more information about our assemblies. We'll run it with a reference this time. 

```bash
# with a ref (hap1)
quast -o 3_quast/ref -t 2 --labels "trio_hap1,trio_hap2,hifi_hap1,hifi_hap2" 2_assembly/trio.hap1.fa 2_assembly/trio.hap2.fa 2_assembly/hifionly.hap1.fa 2_assembly/hifionly.hap2.fa -r 0_data/chr2_35_to42Mbp_hap1.fa
```

❓**Question:**
- Why would it take longer to run QUAST with a reference genome vs without one? And what extra information do you think this might produce?


Open the `html` report in `3_quast/ref` in a web browser and open the QUAST report to answer the following:
- Which assembly has the highest Genome Fraction and what might this mean?
- What does it mean that our assemblies have "misassemblies" when compared with a reference?
- What about mismatches? What would be considered a mismatch? 
- Which assembly is the most similar to the reference? 
- Given that we are comparing assemblies of two separate haplotypes with a single reference sequence, does the higher rate of misassemblies, mismatches and indels in trio haplotype 2 mean that it is of a lower quality than trio haplotype1? Why or why not? 
- Examine the contig length distribution. Why is the contiguity of an assembly not a measure of assembly accuracy but is important to assembly quality? 

Click the "View in Icarus contig browser" and take some time to explore. 
When you're ready, try to answer the following:
- Which haplotype in the hifi assembly is the trio haplotype 1 assembly more similar to?
- Zoom in to see if there are any regions in the other Hifi haplotype that match with the trio haplotype 1 assembly better. If there are, why might this be?

# **Some Final Questions**

These questions will help you better understand some of the concepts covered in the online content that will *probably* help with the upcoming quiz.

- How are haplotype resolved assemblies different from collapsed assemblies? 
- What do you predict would be the differences between an assembly generated from 40x PacBio HiFi reads and one from 40x Nanopore reads? 
- If you had 40x of PacBio HiFi data to assemble a large complex eukaryotic genome, what other types of sequencing data and techniques might you use to produce the best quality assembly possible? 




