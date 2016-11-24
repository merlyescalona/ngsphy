# NGSphy

NGSphy, is a tool designed specifically as an addendum to [`SimPhy`](https://github.com/adamallo/SimPhy)[**Mallo et al. 2015**][^Mallo2015]
a phylogenomic simulator of gene, locus and species trees that considers incomplete lineage sorting, gene duplication and loss and horizontal gene transfer. This tool is able to use SimPhy’s output in order to produce Illumina NGS data from haploid/diploid individuals.

# Getting help

Most common issues, doubts and questions should be solved by reading the [wiki](https://gitlab.com/merlyescalona/ngsphy/wikis/home).

If that is not the case or you find any bug, you can post an issue to this repository,
for reproducibility purposes, with the following files attached:
 - the generated log file
 - the settings file
 - `<simphy_project_name>.command` file


# Obtaining SimPhy NGS wrapper

You can get SimPhy NGS wrapper by cloning the following git repository:

    `git@gitlab.com:merlyescalona/simphy-ngs-wrapper.git`


# Requirements

- Python (2.7 or above)
- Python modules:
    - `argparse`
    - `datetime`
    - `logging`
    - `os`
    - `sys`
    - `numpy`
    - `random`
    - `subprocess`
    - `ConfigParser`
- ART[^art] ([Version ChocolateCherryCake or later](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/))





[^Mallo2015]: Mallo D, de Oliveira Martins L, Posada D (2015) SimPhy: Phylogenomic Simulation of Gene, Locus and Species Trees. Syst. Biol. doi: http://dx.doi.org/10.1093/sysbio/syv082
