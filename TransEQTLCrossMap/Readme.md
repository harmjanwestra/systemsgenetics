# Code to determine cross mapping of genes on other parts of the genome
Trans-eQTLs describe distal regulation of gene expression levels by genetic variants, whereas cis-eQTLs generally describe effects of genetic variants on nearby genes.
Generally, the effect sizes of trans-eQTLs are very weak, because trans-eQTLs are assumed to be indirect effects.
In contrast, being more likely direct effects, cis-eQTLs generally have lower P-values and stronger association statistics. 
However, sometimes associations with suspiciously low P-values are observed in trans-eQTLs. While these could be true effects,
they are more likely some kind of artifact. One of these artifacts can be caused by 'cross-mapping', where parts of the trans-eQTL 
gene is mapping close to the genetic variant. For instance, there might be sequence overlap between the trans-eQTL gene and
an inactive pseudogene copy of the trans-eQTL gene, located in *cis* of the genetic variant. As such, the observed trans-eQTL should
rather be considered a cis-eQTL effect. 

This phenomenon has been described (amongst others) by Saha and Battle: https://pubmed.ncbi.nlm.nih.gov/30613398/.

The code in this repository was used for eQTLgen phase 1 and 2 (https://eqtlgen.org/) and MetaBrain to determine 
whether (parts of) trans-eQTL genes map closely to trans-eQTL variants.

## General procedure
- First we load in a gene annotation GTF file, a list of trans-eQTL variant-gene pairs to test, and a human genome fasta file
- For each eQTL:
  - Either: A) we determine the location of the cis-eQTL gene, and create a reference sequence as follows:
     - Left bound: gene start position - windowsize
     - Right bound: gene stop position + windowsize
     - Window will be clipped if left bound < 0 or right bound > chromosome size
  - Or: A) we determine the location of the variants, and determine the sequence surrounding the variant given a certain window size. We will use this sequence as a genome reference later.
     - The variant will be located in the middle, with some exceptions.
        - Left bound: variant position - windowsize 
        - Right bound: variant position + windowsize
        - Window will be clipped if left bound < 0 or right bound > chromosome size 
  - We next determine the sequence of the trans-eQTL gene, excluding introns, but including all transcripts defined in the annotation.
     - We scan the gene from leftmost to rightmost position, using a certain window size, making sure that each window overlaps the previous window with a certain number of basepairs. 
     - All sequences underneath each window are then written to a fasta file for the gene. 
     - We do this procedure for all possible transcripts defined in the GTF file. 
  - We use BWA to align the generated gene sequence parts to each of the cis-eQTL gene or variant reference sequences
  - We count the aligned basepairs per reference-trans-eQTL gene pair and calculate the proportion of bases aligned over the total number of bases generated for the trans-eGene

The A) procedure for defining a reference window to align trans-eQTL gene reads to is similar to the one presented by Saha and Battle.

The B) procedure is slightly different from the one presented by Saha and Battle, because the procedure above does not assume that the 
gene causing the cross-mapping artifacts is actually annotated as being a gene near the variant. As such, this procedure also can find
potential cross-mapping artifacts for genes that are not supposed to exist. This is somwhat stricter than the gene body approach, and  what we did for eQTLgen and MetaBrain.


## Compilation
This software depends on other modules in the systems genetics repository, mostly from the genetica-libraries module.
Please attempt to compile that module first. The easiest way is to pull the whole systems genetics repository and load them into your favorite IDE as 
a Maven project. Then compile the genetica-libraries module, followed by this module.

## Command line options
| Option | Description |
| ------ | ------------|
|-a,--gtf <arg>|GTF file
|-c,--usecisgeneasrefwindow|Use cis-gene as reference window [default: false]|
|-e,--eqtlfile <arg>|eQTL file|
|-g,--genomefile <arg>| Genome fasta file|
|-i,--indir <arg>|Input directory|
|-m,--mode <arg>|Mode: makeseq or quantify|
|-o,--out <arg>|Output directory or file|
|-r,--readlen <arg>|Read length [default 35]|
|-s,--readshift <arg>|Shift read window with [default 2]|
|-t,--jobtemplate <arg>|Job template file|
|-w,--refwindowsize <arg>|Reference window size [default 5000000]|

## Example of running the software

### Step 1: grepping the sequences from the genome reference

#### Procedure A - reference window is gene body
    java -Xmx16g -jar TransEQTLCrossMap-1.0.0-SNAPSHOT-jar-with-dependencies.jar \
        --mode makeseq \
        --eqtlfile /path/to/eqtlfile.txt \
        --genomefile /path/to/Homo_sapiens.fasta.gz \
        --gtf /path/to/Homo_sapiens.gtf.gz \
        --out /path/to/outputdir/ \
        --refwindowsize 1000 \
        --usecisgeneasrefwindow \
        --readlen 35 \
        --reafshift 10 \
        --jobtemplate jobtemplate.tpl 

#### Procedure B - reference window centered around variant
    java -Xmx16g -jar TransEQTLCrossMap-1.0.0-SNAPSHOT-jar-with-dependencies.jar \
        --mode makeseq \
        --eqtlfile /path/to/eqtlfile.txt \
        --genomefile /path/to/Homo_sapiens.fasta.gz \
        --gtf /path/to/Homo_sapiens.gtf.gz \
        --out /path/to/outputdir/ \
        --refwindowsize 5000000 \
        --readlen 35 \
        --reafshift 10 \
        --jobtemplate jobtemplate.tpl


#### Expected input
**eqtlfile.txt**: A file containing variant/gene pairs with columns for each, tab separated, with header, in the following format:
Expected/required columns:

|Column name| Description                       |
|-----------|-----------------------------------|
|variantchr| chromosome of the variant         |
|variantpos| base pair position of the variant |
|variantname| rsId or other variant name        |
|cisgene| cis gene Ensembl gene ID          |
|transgene| trans gene Ensembl gene Id        |

If using procedure B, the cisgene column is not required and can be replaced with a random string. Make sure that the cis/trans gene IDs used match those defined in the GTF file (e.g. with our without the dot gene version separator).
Note: at this moment, the software does not support X, Y, and MT chromosomes. 

**Homo_sapiens.fasta.gz**: Path to the fasta file containing human genome reference. make sure the build matches the GTF annotation and those of the eQTL results

**Homo_sapiens.gtf.gz**: Path to the GTF file containing gene reference positions. Make sure it contains the locations of the transcripts and exons, and that it matches the build of the reference genome and eQTL input.

**/path/to/output/**: Path to output folder. This is where the alignments and sequences will be stored

**refwindowsize**: Integer value, defines the window size surrounding the genetic variant. When using procedure A, select any reasonable small value (e.g. 1000). For procedure B, use whatever value was used to determine the trans-eQTL window. E..g. for eQTLgen and MetaBrain we used 5000000.

**readlen**: Integer value, size of the sequence parts generated for the gene. We used 35 for eQTLgen.

**reafshift**: Integer value, overlap of sequence between consecutive gene sequence windows. We used 10 for eQTLgen.

**jobtemplate.tpl**: Path to job template file. This is a shell script template that will be used to create jobs for your compute environment. Example provided in this repository.

#### Expected output
In the outdir, several subfolders will be created.

|Subfolder| Description / contents                                                                                         |
|---------|----------------------------------------------------------------------------------------------------------------|
|1-indexjobs| Subfolder containing a number of shell scripts, derived from jobtemplate.tpl. See next step.                   |
|2-alignjobs| Subfolder containing a number of shell scripts, derived from jobtemplate.tpl. See next step.                   |
|3-samsejobs| Subfolder containing a number of shell scripts, derived from jobtemplate.tpl. See next step.                   |
|alignments| Folder where eventual alignments will be stored                                                                |
|refwindows| Folder where .fa.gz files are stored for the reference sequences for each gene-body or variant-centered window |
|testgenes| Folder where .fa.gz files are stored with 'read' sequences generated from the trans-eQTL gene.                 |

### Step 2: indexing the reference sequences and alignment
Step 1 created the folders 1-indexjobs, 2-alignjobs, and 3-samsejobs.
First, run the shell scripts in the 1-indexjobs folder, then those in the 2-alignjobs folder, and finally, those in the 3-samsejobs folder.
These jobs can also be submitted to an HPC environment. Each shell script contains a number of BWA commands. Total runtime should be about 
15 minutes total, but depends on the number of indexing and the number of alignments to be performed.
BWA alignments should ultimately appear in the alignments subfolder, recognizable by the .sai and .sam file extensions. 

### Step 3: quantify cross-mapping
In this final step, we'll use the output generated above to quantify the amount of potential cross mapping.

#### Procedure A - reference window is gene body

    java -Xmx16g -jar TransEQTLCrossMap-1.0.0-SNAPSHOT-jar-with-dependencies.jar \
        --mode quantify \
        --eqtlfile /path/to/eqtlfile.txt \
        --gtf /path/to/Homo_sapiens.gtf.gz \
        --indir /path/to/outdirofstep1/ \
        --out /path/to/outputfile.txt \
        --refwindowsize 1000 \
        --usecisgeneasrefwindow

#### Procedure B - reference window centered around variant
    java -Xmx16g -jar TransEQTLCrossMap-1.0.0-SNAPSHOT-jar-with-dependencies.jar \
        --mode quantify \
        --eqtlfile /path/to/eqtlfile.txt \
        --gtf /path/to/Homo_sapiens.gtf.gz \
        --indir /path/to/outdirofstep1/ \
        --out /path/to/outputfile.txt \
        --refwindowsize 5000000

