# SimPhy NGS wrapper

This is a plugin for [SimPhy](https://github.com/adamallo/SimPhy) [ **Mallo et al. 2015**][^Mallo2015] *A comprehensive simulator of gene family
evolution*. SimPhy NGS Wrapper generates diploid individuals from the sequences
generated of a SimPhy project, and afterwards generates reads from such
individuals with a next-generation sequencing simulator, ART.

For more information about usage and installation please go to the [wiki page](https://gitlab.com/merlyescalona/simphy-ngs-wrapper/wikis/home)

## General:

It accepts multiple prefixes (SimPhy data set prefixes) to do the mating. At
least one must be given. All the datasets will be mated in the same manner. The
mating information is written down inside the <SimPhy_project> folder in a file
 under the name of:

    <SimPhy_project>.mating

The file is in CSV format. Columns of the file represent:

    indexST,indexLOC,individualID,speciesID,mateID1,mateID2

When iterating over the species tree replicates, if a replicate does not have
the corresponding sequences, such replicate will be skipped.

All the FASTA files generated are going to be stored in the SimPhy project
folder under the "individuals" folder. There are multiple options to generate
the individual files. It is possible to get 1, 2 or 3 files per individual.
When selecting the number of files to output:

  1. Only one file for both sequences.
  2. A file per strand equence.
  3. A file with both sequences and a file per strand sequence.

By default, the option 3 will be used.

Name convention for this files is the following:

- Project_ST_GT_PREFIX_IND.fasta
- Project_ST_GT_PREFIX_IND_DES_R1.fasta
- Project_ST_GT_PREFIX_IND_DES_R2.fasta

Where:
- Project: SimPhy's project name
- ST: Index of the species tree used to generate the specific individual.
- GT: Index of the gene tree used to generate the specific individual.
- PREFIX: Prefix to which corresponds the sequence used for the individual.
- IND: ID of the individual created for a specific gene tree.
- DES:  When a file is generated per strand of an individual, here is the
        description of the sequence, this description specifies one of
        the sequences used to generate the individual.
- S1/S2:    Tag that indetifies the strand of the individual.


[^Mallo2015]: Mallo D, de Oliveira Martins L, Posada D (2015) SimPhy: Phylogenomic Simulation of Gene, Locus and Species Trees. Syst. Biol. doi: http://dx.doi.org/10.1093/sysbio/syv082
