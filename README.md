# SimPhy NGS wrapper

## Description:

Program to simulate mating. Extracts sequences from SimPhy project folder and
trees to generate diploid individuals.

## Assumptions:

- There is always one (1) outgroup. The outgroups will be considered to have a
single sequence, though the outputted individual will have such sequence
duplicated.

- The number of gene trees is the same for all species tree replicates.
- The number of datasets should be the same for all species tree replicates.

## Installation:

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
