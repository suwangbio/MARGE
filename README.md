# MARGE
### Model-based Analysis of Regulation of Gene Expression

## Summary
MARGE is a robust methodology that leverages a comprehensive library of genome-wide H3K27ac ChIP-seq profiles to predict key regulated genes and cis-regulatory regions in human or mouse. The framework has three main functions: MARGE-potential, MARGE-express, and MARGE-cistrome.

## Introdunction
#### MARGE-potential 
MARGE-potential defines a regulatory potential (RP) for each gene as the sum of H3K27ac ChIP-seq signals weighted by a function of genomic distance from the transcription start site. The MARGE framework includes a compendium of RPs derived from 365 human and 267 mouse H3K27ac ChIP-seq datasets. Relative RPs, scaled using this compendium, are superior to super-enhancers in predicting BET-inhibitor repressed genes.

#### MARGE-express 

MARGE-express, uses regression to link gene expression perturbations with regulatory potentials derived from a small subset of H3K27ac ChIP-seq data from the full compendium. In this way MARGE determines changes in regulatory potentials that are predictive of gene expression changes. Many biological perturbations that have been profiled using microarray or RNA-seq expression analysis have not been studied using chromatin profiling techniques such as ChIP-seq or DNase-seq. MARGE-express serves to identify relationships between gene sets and cis-regulatory environments to enable the study of the cis-regulation of these gene sets. It can even predict differentially regulated genes, such as ncRNAs, that are absent from the platform used to measure gene expression.

#### MARGE-cistrome 

MARGE-cistrome, to predict co-regulated sets of cis-elements. These cis-elements are predicted on the basis of 1kb H3K27ac ChIP-seq windows centered on the union of DNase-seq peaks that have been identified in a wide spectrum of cell types. MARGE-cistrome identifies patterns of perturbations in the cis-elements that are consistent with perturbations in the H3K27ac regulatory potentials identified by MARGE-express. Investigators who wish to understand how particular genes in their gene set are regulated can use MARGE-cistrome to identify candidate cisregulatory elements, even when chromatin profiling data is not available in their system.

## How to install?
Before running MARGE, make sure the following software is installed in system:

1. Python2.7 or newer; Or Python3.4 or newer
2. snakemake (>=3.4.1 is strongly recommended)
3. HDF5 ( >= 1.8.7 is strong recommended)
4. MACS2 ( >= 2.1.0 is strong recommended)
5. Some UCSC tools

bedClip, bigWigToBedGraph, bigWigSummary, and bigWigAverageOverBed 
Download the above four programs and simply put them in a directory (you might need to change the path for these programs in the configure file described below)

#### python dependencies
1. Install setuptools and argparse
2. Install numpy
3. Install scipy
4. Install sklearn
5. Install tables
6. Install twobitreader

#### Install MARGE

```
$ unzip Py2_MARGE.zip (Py3_MARGE.zip)
$ cd Py2_MARGE (Py3_MARGE)
$ python setup.py install --prefix=/path/to/where/to/install
```

## How to use?
### Input files
#### Bam file
If you have a H3K27ac ChIP-seq dataset, and want to know the key genes in this cell or tissue or some specific conditions, you might want to try MARGE to help you know much about the key regulated genes (or master regulators)

When you firstly got the H3K27ac ChIP-seq raw data, it is in FASTQ format, we recommend you use BWA (version 0.7.10) do the mapping, you can also use other aligners like BOWTIE.

Once you have mapped your FASTQ file to the reference genome, you will normally end up with a SAM or BAM alignment file. SAM stands for Sequence Alignment/Map format, and BAM is the binary version of a SAM file. We recommend users using BAM format alignment as the input of MARGE, you can use samtools to convert SAM format to BAM format


1. **Make sure the BAM files with '.bam' as the suffix.** 
2. **If you want to process multiple bam files at the same time, put them in the same directory.**

If bam file is hg38 or mm10, MARGE will output **Regulatory Potential, Relative Regulatory Potential**

If bam file is hg19 or mm9, MARGE will output **Regulatory Potential**
#### Gene List file
Suppose you have a list of interested genes, and want to detect the cis-regulatory regions which might potentially regulate these genes expression. In such case, you can run MARGE with the following input format

MARGE support two kinds of Gene List format: Gene_Only and Gene_Response.

1. **Both Gene_Only and Gene_Response format files should with '.txt' as the suffix.**
2. **If you want to process multiple gene list files at the same time, put them in the same directory.**

**Gene_Only** format, one column, each row is a gene, both GeneSymbol ID and RefSeq ID are supported by MARGE

###### Example
|geneID|
|---|
|AR|
|BAGE4|
|CST1|
|DLX3|
|...|

**Gene_Response** format, two columns, each row is a gene, first column is gene ID and second column is 1(target) or 0(non-target). Both GeneSymbol ID and RefSeq ID are supported by MARGE
###### Example
| geneID |1/0|
|---|---:|
| NM_000397 |1|
| NR_045675 |0|
| NM_033518 |1|
| NM_052939 |1|
| NM_002290 |0|
| NM_000495 |0|
|...|...|
**Please do not contain the header in the gene list file, and make sure it is tab delimitated in the Gene_Response file format.**

**!! Only when assembly is hg38 or mm10, MARGE can do cis-regions prediction.**

**Output: The 10 most relevant H3K27ac samples, which can best interpret these genes expression status, then MARGE use these 10 samples information to predict cis-regions that might potentially regulate the gene expression status**

## Run MARGE
All steps below have to be executed in a terminal.


#### Testing MARGE

We provide some test data and config.json file on the Download page

For usage of MARGE on real data, please see all steps below.

#### Step 1: Choosing a workflow directory

To predict the key regulated genes or cis-regulatory regions in human or mouse with MARGE. First you should select a directory where the workflow shall be executed and results will be stored. Please choose a meaningful name and ensure that the directory is empty.

#### Step 2: Initializing a new workflow

We assume that you have chosen a workflow directory (here path/to/my/workflow), and that you have a set of BAM files (sample1.bam sample2.bam sample3.bam sample4.bam...) in one directory (path/to/marge_testdata/) or a set of Gene List files (genelist1.txt genelist2.txt genelist3.txt...) in one directory (path/to/marge_testdata/) or both of them you want to analyze with MARGE. You can initialize the workflow with

`$ marge init path/to/my/workflow`

This will initialize the MARGE workflow in a given directory. It will install a Snakefile, and a config file in this directory. Configure the config file according to your needs.

#### Step 3: Configure the workflow

Now change the parameters and paths in the config.json manually via

`$ cd path/to/my/workflow`

and open the config.json file. Edit the config file to your needs. Especially, define ASSEMBLY for this job, define MARGEdir, which is the path for the MARGE source code directory, define REFdir for the path to reference directory which you can download from MARGE Library, SAMPLESDIR and SAMPLES for the BAM samples, EXPSDIR and EXPS for the Gene List samples, and EXPTYPE for the format of Gene List file and ID type (RefSeq or GeneSymbol) for the gene used in Gene List. Change the path to some tools like MACS2, bedGraphToBigWig, bedClip, bigWigAverageOverBed, and bigWigSummary 

#### Step 4: Execute the workflow

Once configured, the workflow can be executed with Snakemake. First, it is advisable to invoke a dry-run of the workflow with

`$ snakemake -n `

which will display all jobs that will be executed. If there are errors in your config file, you will be notified by Snakemake. The actual execution of the workflow can be started with

`$ snakemake --cores 8`

here using 8 CPU cores in parallel. For more information on how to use Snakemake, refer to the Snakemake Documentation.
