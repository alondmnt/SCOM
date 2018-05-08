## SCOM algorithm

This repository mirrors the original paper's [web site](http://www.cs.tau.ac.il/~tamirtul/SCOMs/).

SCOMs - or Spatially Co-evolving Orthologous Modules - are a graph theory-based model for studying the organization of genes in multiple-organisms. These are sub-graphs that comprise of gene families and describes the relations between them . This method was originally used to study the evolution of gene organization (more specifically, their 3D distribution inside the nucleus*) in the budding yeast _S. cerevisiae_ and the fission yeast _S. pombe_. Moreover, it can easily be extended for additional purposes, such as comparing different tissues from the same organism, monitoring the same cells under various conditions (that may affect the 3D organization of chromosomes), and for a differential study of healthy and cancer cells.

\* using chromosome conformation capture (3C, Hi-C) data.

### Cite

Diament A, Tuller T. Tracking the evolution of 3D gene organization demonstrates its connection to phenotypic divergence. Nucleic Acids Research, 2017. [DOI:10.1093/nar/gkx205](http://dx.doi.org/10.1093/nar/gkx205)

### SCOM scripts

The scripts here require MATLAB (for generating conservation / divergence networks from Hi-C data) and python (for detecting SCOMs in any given network(s), tested on python 3). The provided shell script 'example.sh' demonstrates the pipeline for a small-scale example.

### SCOM database

For an online interactive database of annotated conserved and divergent SCOMs between the budding yeast _S. cerevisiae_ and the fission yeast _S. pombe_, see the [NDEx site](http://www.ndexbio.org/#/user/177f1419-f292-11e6-a7f1-0ac135e8bacf).

