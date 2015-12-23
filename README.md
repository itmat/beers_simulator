

# beers_simulator

Benchmarker for Evaluating the Effectiveness of RNA-Seq Software (**BEERS**)

**BEERS** is a simulation engine for generating RNA-Seq data.

**BEERS** was designed to benchmark RNA-Seq alignment algorithms and also algorithms that aim to reconstruct different isoforms and alternate splicing from RNA-Seq data.

**Publication:** [Comparative Analysis of RNA-Seq Alignment Algorithms and the RNA-Seq Unified Mapper (RUM)](http://www.ncbi.nlm.nih.gov/pubmed/21775302?dopt=Abstract)
**Authors:** Gregory R. Grant, Michael H. Farkas, Angel Pizarro, Nicholas Lahens, Jonathan Schug, Brian Brunk, Christian J. Stoeckert Jr, John B. Hogenesch and Eric A. Pierce.

By default **BEERS** simulates either mouse or human paired-end RNA-Seq data modeled on the illumina platform. It starts with a large number of gene models (approx 500K) taken from about ten different published annotation efforts, and then chooses a fixed number of these genes at random (30,000 by default). This avoids biasing for or against any particular set of annotations. **BEERS** then introduces substitutions, indels, alternate spice forms, sequencing errors, and intron signal. **BEERS** can also simulate strand specific reads. **BEERS** does not simulate quality scores. There are four configuration files required (available below).

**BEERS** can also be configured to use any set of gene models. Pre-built indexes for human refseq are given below. Using these indexes will generate a much tamer set of transcripts.

**BEERS** is written in perl. See below for more information and run with no parameters to get the usage. *Math::Random* must be installed on your system.

## Configuration files

**BEERS** requires four configuration files, these are available for human and mouse:

 - [Mouse config files](http://itmat.rum.s3.amazonaws.com/simulator_config_mouse.tar.gz) The files are based on eleven annotation tracks: UCSC, RefSeq, RefSeq-Other, Ensembl, Vega, AceView, GenScan, GeneID, NSCAN, SGP and Transcriptome.
 -  [Human config files](http://itmat.rum.s3.amazonaws.com/simulator_config_human.tar.gz) The files are based on ten annotation tracks: UCSC, RefSeq, RefSeq-Other, Ensembl, Vega, AceView, GenScan, GeneID, NSCAN, and SGP
 - [Human RefSeq config files](http://itmat.rum.s3.amazonaws.com/simulator_config_refseq.tar.gz) This one can be used to create a much tamer set of transcripts, that relies solely on refseq annotation.
To use these refseq specific config files, run with the "-configstem refseq" option, and also set "-palt 0". This will generate simulated data based solely on refseq, without introducing further alternate forms than what already exist in refseq (that's what the -palt 0 option does). By default intron signal will also be generated, you can turn that off with the -nointron option.


<img align="right" src="http://www.cbil.upenn.edu/BEERS/simulator_workflow_small.jpg">
In order not to bias for or against any particular set of gene models, many different sets of annotation were merged (e.g. AceView, Ensembl, Geneid, Genscan, NSCAN, RefSeq, SGP, Transcriptome, UCSC, Vega). These models were filtered to remove most of the junctions that had uncharacterized splice signals, The simulator workflow is shown below. The first module chooses N of the gene models at random, with a default of N=30,000. This is done in order to not bias our data towards any particular set of gene models. Alternate splice forms are then created for each gene by preferentially leaving in exons, where the number of alternate forms is a parameter with a default of two. The percentage of signal coming from alternate splice forms is a parameter with a default of 20%. Polymorphisms (indels and substitutions) are introduced into the exons, according to independent rates. A gene quantification file is used to assign an empirical distribution of signal that mimics real data. This file is further used to determine the distribution of intronic signal, so that preferential intron inclusion can be simulated. Reads are then produced by choosing a gene at random, possibly leaving in an intron, choosing a fragment of normally distributed length, introducing random base and tail error, and then reporting the M bases of the fragment from either end, where M is the read length. Random base error is set according to one parameter and tail error is set according to three parameters: the percent of low quality tails, the length of the low quality tail, and the quality of the low quality tail. The reads generated are reported to a fasta file. The true coordinates of each real and the true junctions spanned are reported to text files. The set of gene models used, the alternate splice forms, and the polymorphisms are reported to log files.






## Descriptions of the output files

 1. file ending .cig has the true alignments in several formats for easy parsing
 2. file ending in .fa has the actual reads, forward are the "a" reads, reverse the "b" reads. 10,000,000 pairs.
 3. simulated_reads_junctions-crossed.txt
   - this file has the information for each read which crossed a junction
   - you need this because when calculating the false negative rate, it's not fair to any algorithm to use junctions that weren't actually crossed by any reads.
 4. simulated_reads_transcripts.txt
   - this file has the full transcript info for each simulated transcript (basically the same as refseq with some redundancy issues fixed and zero length introns removed).
   - you need this to calculate the false positive rates.

There are a few other files, basically all the info about the 'truth' you could possibly want to know is in there somewhere including the coverage plot and the fasta file of genes and the signal intensities.
