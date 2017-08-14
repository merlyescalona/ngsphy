NGSphy: simulation of genome-wide next-generation sequencing data from species trees
================================================================================

About NGSphy
--------------------------------------------------------------------------------

NGSphy, is a tool designed as a pipeline to generate next-generation sequencing
(NGS) data from genome sequences of haploid/diploid individuals underlain by gene
trees or species trees. This tool was originally thought as an addendum to: SimPhy
(https://github.com/adamallo/simphy) a phylogenomic simulator of gene, locus and
species trees that considers incomplete lineage sorting, gene duplication and
loss and horizontal gene transfer, to be able to handle its output and generate
data from it. NGSphy is able to use SimPhy’s output or INDELible's input in
order to generate haploid or diploid individuals and/or produce Illumina data
using, ART[^art], or read counts.

What can be done with NGSphy
--------------------------------

- Generate haploid individuals from species tree distributions
- Generate diploid individuals from species tree distributions
- Generate genome sequences of haploid individuals from a random reference sequence and a given gene tree
- Generate genome sequences of diploid individuals from a random reference sequence and a given gene tree
- Generate genome sequences of haploid individuals from a given reference sequence and a given gene tree
- Generate genome sequences of diploid individuals from a given reference sequence and a given gene tree
- Generate NGS Illumina reads of genome sequences of haploid individuals
- Generate NGS Illumina reads of genome sequences of diploid individuals
- Generate NGS read counts of genome sequences of haploid individuals
- Generate NGS read counts of genome sequences of diploid individuals
- For the NGS data generation, variation of coverage due to the following:
    - genomic stochasticity - noise
    - variation across individuals and/or loci
    - phylogenetic decay
    - on/off target effect from targeted-sequencing experiments.

Citation
===============


If you use NGSphy, please cite:

- Escalona, M., Rocha S. and Posada D. (XXXX) *NGSphy: simulation of genome-wide next-generation sequencing data from species trees*

- Weichun Huang, Leping Li, Jason R Myers, and Gabor T Marth (2012) *ART: a next-generation sequencing read simulator*. Bioinformatics  28 (4): 593-594

- William Fletcher and Ziheng Yang  (2009) *INDELible: A flexible simulator of biological sequence evolution*. Molecular Biology and Evolution. 26 (8): 1879–88. [doi:10.1093/molbev/msp098](doi:10.1093/molbev/msp098).

- Diego Mallo, Leonardo De Oliveira Martins and David Posada (2015). *SimPhy : Phylogenomic Simulation of Gene, Locus, and Species Trees*. Systematic Biology., November, syv082. [doi:10.1093/sysbio/syv082](doi:10.1093/sysbio/syv082).

- Sukumaran, J. and Mark T. Holder. 2010. *DendroPy: A Python library for phylogenetic computing*. Bioinformatics 26: 1569-1571.


Obtaining NGSphy and dependencies
================================================================================

You can get NGSphy by cloning the following git repository:

- git@gitlab.com:merlyescalona/NGSphy.git

Third-party apps:

- ART: Version ChocolateCherryCake or later(http://www.niehs.nih.gov/research/resources/software/biostatistics/art/)
- INDELible:  Version 1.0.3 (http://abacus.gene.ucl.ac.uk/software/indelible/)
- INDELible:  Version  (own/github)

- Python (2.7) (It has only been tested with this version)
- Python modules:
  - argparse
  - ConfigParser
  - datetime
  - logging
  - numpy
  - os
  - random
  - scipy
  - subprocess
  - sys


For more information please go to the wiki.

Getting help
===============
Most common issues, doubts and questions should be solved by reading the wiki (https://gitlab.com/merlyescalona/ngsphy/wikis/home).

If that is not the case or you find any bug, you can post an issue to this repository,
for reproducibility purposes, with the following files attached:
 - the generated log file
 - the settings file
 - <simphy_project_name>.command file

Usage
===============

::
usage: ngsphy.py  [-s <settings_file_path>]
               [-l <log_level>] [-v] [-h]

- **Optional arguments:**
    - -s <settings_file_path>, --settings <settings_file_path>: Path to the [settings](./2-Settings) file.
    - -l <log_level>, --log <log_level>: Specified [log levels](./2-Settings#log-levels) that will be shown through the standard output. Entire log will be stored in a separate file.
       - Values:[DEBUG, INFO, WARNING,ERROR].
       - Default: INFO.

- **Information arguments:**
    - -v, --version: Show program's version number and exit
    - -h, --help: Show help message and exit

Third-party software involved
================================================================================

SimPhy
--------------------------------------------------------------------------------

SimPhy is a program for the simulation of gene family evolution under incomplete lineage sorting (ILS), gene duplication and loss (GDL), replacing horizontal gene transfer (HGT) and gene conversion (GC). SimPhy simulates species, locus and gene trees with different levels of rate heterogeneity, and uses INDELible to evolve nucleotide/codon/aminoacid sequences along the gene trees. The input for SimPhy are the simulation parameter values, which can be fixed or sampled from user-defined statistical distributions. The output consists of sequence alignments and a relational database that facilitate posterior analyses.

INDELible
--------------------------------------------------------------------------------

INDELible is an application for biological sequence simulation that combines many features Using a length-dependent model of indel formation it can simulate evolution of multi-partitioned nucleotide, amino-acid, or codon data sets through the processes of insertion, deletion, and substitution in continuous time.

Nucleotide simulations may use the general unrestricted model or the general time reversible model and its derivatives, and amino-acid simulations can be conducted using fifteen different empirical rate matrices. Substitution rate heterogeneity can be modelled via the continuous and discrete gamma distributions, with or without a proportion of invariant sites. INDELible can also simulate under non-homogenous and non-stationary conditions where evolutionary models are permitted to change across a phylogeny.

Unique among indel simulation programs, INDELible offers the ability to simulate using codon models that exhibit nonsynonymous/synonymous rate ratio heterogeneity among sites and/or lineages.


ART
--------------------------------------------------------------------------------

ART is a set of simulation tools to generate synthetic next-generation sequencing reads. ART simulates sequencing reads by mimicking real sequencing process with empirical error models or quality profiles summarized from large recalibrated sequencing data. ART can also simulate reads using user own read error model or quality profiles. ART supports simulation of single-end, paired-end/mate-pair reads of three major commercial next-generation sequencing platforms: Illumina’s Solexa, Roche’s 454 and Applied Biosystems’ SOLiD. ART can be used to test or benchmark a variety of method or tools for next-generation sequencing data analysis, including read alignment, *de novo* assembly, SNP and structure variation discovery. ART outputs reads in the FASTQ format, and alignments in the ALN format. ART can also generate alignments in the SAM alignment or UCSC BED file format.

Motivation
================================================================================
This has been made to be included in a pipeline for the simulation of NGS-like data from phylogenies of closely related species, based on a targeted sequencing experiment.
The pipeline includes:

1. Simulation of species and gene trees using SimPhy.
2. Simulation of captured loci using INDELible[^indelible]. (varying size of the loci according to the probes designed)
3. Generation of the individuals (haploid/diploid)
4. Simulation of NGS data from the independent loci, whether Illumina reads or read counts.


# Data origin
================================================================================

NGSphy generates sequencing data from haploid or diploid individuals within species trees and/or gene trees. Data origin can be selected from one of the following:

1. species tree distributions (a SimPhy project)
2. a Newick file with an INDELible control file (INDELible **Input A**)
3. a reference sequence, a Newick file and an INDELible control file (INDELible **Input B**)

SimPhy project
--------------------------------------------------------------------------------

A valid SimPhy project has been obtained from a complete SimPhy run, which consists on the execution of SimPhy to simulate species
trees and gene trees; the generation of the INDELible control files according to the parameters established (running of the INDELible_wrapper.pl
script) and the posterior genome sequence simulation, through INDELible. Detailed description of SimPhy's output can be found
[here](https://github.com/adamallo/SimPhy/wiki/Manual#53-output-files)). In general terms, the  output is structured as:

- A main folder stores all the results.
    - a sub-folder per replicate and a set of files with information of the simulations, within each replicate folder we will find basically
     the species, locus and gene trees as well as the genome sequences for each gene tree.
    - Related to the files with information of the simulations, for running NGSphy we are interested in only three (3) of the files:
        - A plain text file with the original command line arguments (X.command)
        - A plain text file summarizing the sampled options (X.params).
        - A SQLite database composed by three (3) linked tables with different information about the species, locus and gene trees.

INDELible A : random reference
--------------------------------------------------------------------------------

This will represent the generation of genome sequences from a random reference sequence. This INDELible dataset (**Input A**) consist on a
given gene tree in Newick format and a control file with the parameters of the corresponding evolution model.

INDELible B : specific reference
--------------------------------------------------------------------------------

This will represent the generation of genome sequences from a given reference sequence. For this data origin type, INDELible**Input B**,
it is necessary to have a gene tree in Newick format, a control file with the parameters of the corresponding evolution model, a reference sequence and
select which tip, from the given gene tree, will be the new root of the tree.


Process
================================================================================

NGSphy, verifies all the content of the project, the settings files involved and/or the existence of the corresponding third-party applications
in order to run. More information about each specific ***Data origin***, [here](./3-data-origin).

If the data origin corresponds to the INDELible configuration with a given reference sequence, first thing that will be done is a re-rooting of
the given tree, to the selected tip. Afterwards, independently of the given or random sequence, INDELible will be ran with the given configuration
parameters to obtain the expected genome sequences.

Then, there will be the *generation of individuals*, whether haploid or diploid:

- For haploid individuals, resulting genome sequences are separated into single FASTA files and identified. In addition, a file is generated with the
relationship existing within the individual generated and the description of the sequence it belongs to.

- For diploid individuals, there is a process of verification that the project content includes species trees with an even number
of individuals per taxa, those species trees can be filtered. Take into account that it will make no sense to generate diploid
individuals and leave sequences unpaired. Sequences are then "mated", and individuals are generated by randomly sampling without
replacement two sequences within the same gene family. Output will include a mating table for each species tree replicate
(case of SimPhy data) or output folder (INDELible data),with the identifiers for the sequences mated and the individuals generated.

Afterwards, the coverage distribution matrices will be computed according to the parameters introduced (detailed process [here](./6-coverage)).
And finally, the sequencing data generation, consist on either Illumina reads ([ngs-reads-art] section) or read counts (VCF files, [read-count] section).

- For the Illumina reads, program calls out ART, the NGS simulator, with the parameters established in the settings file and generates
reads from the previously generated individuals. Resulting files depend on the settings introduced, they are related to the execution of
the ART processes (scripts and text files) and the output of such processes, basically, ALN,BAM and/or FASTQ files.

- For read counts, two scenarios are involved, the TRUE read counts (no sequencing error) and the SAMPLED (where sequencing error is introduced).
Process consist first on the identification of the variable sites, calculation of the coverage per variable position, which is sampled from a
Negative Binomial distribution (with mean c_{i,l} and equal dispersion) and in the case of diploid individuals, coverage is distributed (again)
per strand following a binomial distribution (with 0.5 probability per strand). Afterwards, haplotypes/genotypes are computed, as well as the
corresponding likelihoods. The output consist then on a set of VCF files, one per locus.
