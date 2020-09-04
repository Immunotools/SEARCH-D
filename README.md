# SEARCH-D
An algorithm for de novo finding immunoglobulin diversity (D) genes in IGH locus. 

# Availability
The tool is available as Supplemental Code of the manuscript "V(DD)J recombination is an important and evolutionary conserved mechanism for generating antibodies with unusually long CDR3s" by Safonova and Pevzner, 2020 at the [Genome Research](https://genome.cshlp.org/) (link will be added upon the publication). 

## Requirements:
- python 2 or 3
- Biopython
- NumPy
- pandas

## Run:
python search_d.py IGHD_locus.fasta output_fname.fasta

## Examples:
* search for IGHD genes in the human IGHD locus
```
python search_d.py IGHD_loci/homo_sapiens_IGHD.fasta human_d_denovo.fasta 
```
* search for IGHD genes in the mouse IGHD locus
```
python search_d.py IGHD_loci/mus_musculus_igh_locus.fa mouse_d_denovo.fasta 
```
* search for IGHD genes in the rat IGHD locus
```
python search_d.py IGHD_loci/RatNor_igh_locus.fasta rat_d_denovo.fasta
```
* search for IGHD genes in the cow IGHD locus
```
python search_d.py IGHD_loci/BosTaurus_igh_locus.fasta cow_d_denovo.fasta
```
