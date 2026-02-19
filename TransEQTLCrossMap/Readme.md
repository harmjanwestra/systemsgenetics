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

The code in this repository was used for eQTLgen phase 1 and 2 (https://eqtlgen.org/) to determine 
whether (parts of) trans-eQTL genes map closely to trans-eQTL variants.

## General procedure
- First we load in a gene annotation GTF file, a list of trans-eQTL variant-gene pairs to test, and a human genome fasta file
- For each variant-gene pair:
  - We determine the location of the variants, and determine the sequence surrounding the variant given a certain window size. We will use this sequence as a genome reference later.
    - The variant will be located in the middle, with some exceptions.
      - Left bound: (variant position - (windowsize/2)) 
      - Right bound  (variant position + (windowsize/2))
      - Window will be clipped if left bound < 0 or right bound > chromosome size
  - We next determine the sequence of the gene, while choosing to include or exclude exons.
    - We scan the gene from leftmost to rightmost position, using a certain window size, making sure that each window overlaps the previous window with a certain number of basepairs. 
    - All sequences underneath each window are then written to a fasta file for the gene. 
    - We do this procedure for all possible transcripts defined in the GTF file. 
  - We use BWA to align the generated gene sequence parts to each of the variant reference sequences
  - We count the aligned basepairs per variantreference-gene pair and calculate the proportion of bases aligned over the total number of bases generated for the gene

This procedure is slightly different from the one presented by Saha and Battle, because the procedure above does not assume that the 
gene causing the cross-mapping artifacts is actually annotated as being a gene near the variant. As such, this procedure also can find
potential cross-mapping artifacts for genes that are not supposed to exist.


## Compilation
This software depends on other modules in the systems genetics repository, mostly from the genetica-libraries module.
Please attempt to compile that module first. The easiest way is to pull the whole systems genetics repository and load them into your favorite IDE as 
a Maven project. Then compile the genetica-libraries module, followed by this module.

## Running the software

### Step 1: grepping the sequences from the genome reference
    java -Xmx16g -jar TransEQTLCrossMap-0.1.0-SNAPSHOT-jar-with-dependencies.jar \
        makeseq \
        /path/to/transeqtls.txt \
        /path/to/Homo_sapiens.fasta.gz \
        /path/to/Homo_sapiens.gtf.gz \
        /path/to/output/ \
        ciswindowsize[int] \
        readwindowsize[int] \
        readwindowoverlap[int] \
        exportgenes[bool] \
        exportexons[bool] \
        jobtemplate.tpl \
        [mode:bulk|5prime|3prime]

#### Expected input
**transeqtls.txt**: A file containing variant/gene pairs with columns for each, tab separated, with header, in the following format:
EnsemblGeneID\tSNPId\tSNPChr\tSNPPosition. For example:

    GeneID   VariantId  VariantChr  VariantPos
    ENSG000001  rs123  1  1001

Note: at this moment, the software does not support X, Y, and MT chromosomes. Make sure that the gene IDs match the ones defined in the GTF file

**Homo_sapiens.fasta.gz**: Path to the fasta file containing human genome reference. make sure the build matches the GTF annotation and those of the eQTL results

**Homo_sapiens.gtf.gz**: Path to the GTF file containing gene reference positions. Make sure it contains the locations of the transcripts and exons, and that it matches the build of the reference genome and eQTL input.

**/path/to/output/**: Path to output folder. This is where the alignments and sequences will be stored

**ciswindowsize**: Integer value, defines the window size surrounding the genetic variant. We used 10000000 for eQTLgen

**readwindowsize**: Integer value, size of the sequence parts generated for the gene. We used 35 for eQTLgen

**readwindowoverlap**: Integer value, overlap of sequence between consecutive gene sequence windows. We used 10 for eQTLgen

**exportgenes[bool]**: Boolean (true/false) value: export each gene as separate fasta file. We used true for eQTLgen.

**exportexons[bool]**: Boolean (true/false) value: export each exon as a separate fasta file. We used false for eQTLgen.

**jobtemplate.tpl**: Path to job template file. This is a shell script template that will be used to create jobs for your compute environment. Example provided in this repository.

**[mode]**: Optional argument; possible settings: bulk, 5prime, 3prime. Choose between generating sequences for the full gene (bulk), only the 3' end, or only the 5' end. The latter two options were created to make the generated sequences more similar to those observed in a single cell experiment.  
