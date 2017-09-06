===============
NGSphy
===============

© 2017 Merly Escalona (<merlyescalona@uvigo.es>), Sara Rocha, David Posada

University of Vigo, Spain ([http://darwin.uvigo.es](http://darwin.uvigo.es))

About NGSphy
===============
NGSphy is a Python open-source tool for the genome-wide simulation of NGS data (read counts or Illumina reads) obtained from thousands of gene families evolving under a common species tree, with multiple haploid and/or diploid individuals per species, where sequencing coverage (depth) heterogeneity can vary among species, individuals and loci, including off-target and untargeted loci.

Documentation
===============

The documentation is available  here (https://github.com/merlyescalona/ngsphy/wiki)

Citation
===============

If you use NGSphy, please cite:

- Escalona, M, Rocha S and Posada D. *NGSphy: phylogenomic simulation of NGS data*. Submitted.

- if running ART cite also:
    - Huang W, Li L, Myers JR and Marth, GT. (2012) *ART: a next-generation sequencing read simulator*. Bioinformatics  28 (4): 593-594

- if using SimPhy cite also:
    - Mallo D, De Oliveira Martins L and Posada D. (2016). *SimPhy : Phylogenomic Simulation of Gene, Locus, and Species Trees*. Systematic Biology 65(2): 334-344.

- if using single gene tree inputs, cite also:
    - Fletcher, W and Yang Z. (2009) *INDELible: A flexible simulator of biological sequence evolution*. Molecular Biology and Evolution. 26 (8): 1879–88.
    - Sukumaran, J and Holder MT. (2010). *DendroPy: A Python library for phylogenetic computing*. Bioinformatics 26: 1569-1571.

Input/output files
===============

Input
---------------

[Single gene-tree scenario]

- NGSPhy settings file
- INDELible control file
- Newick file with single gene tree
- ancestral sequence file (optional)
- reference allele file (optional)

[Species-tree scenario]

- NGSPhy settings file
- SimPhy output
- reference allele file (optional)

Output files
---------------
- NGS reads:
    - FASTQ
    - ALN
    - BAM
- read counts:
    - VCF
- sequence alignments:
    - FASTA
- coverage variation
    - CSV
- log files
- bash scripts
