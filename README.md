# SEARCH-D
An algorithm for de novo finding IGHD genes in IGH locus

Requirements:
- python 2 or 3
- Biopython
- NumPy
- pandas

Run:
python search_d.py IGHD_locus.fasta output_fname.fasta

Examples:
python search_d.py IGHD_loci/homo_sapiens_IGHD.fasta human_d_denovo.fasta # search for IGHD genes in the human IGHD locus
python search_d.py IGHD_loci/mus_musculus_igh_locus.fa mouse_d_denovo.fasta # search for IGHD genes in the mouse IGHD locus
python search_d.py IGHD_loci/RatNor_igh_locus.fasta rat_d_denovo.fasta # search for IGHD genes in the rat IGHD locus
python search_d.py IGHD_loci/BosTaurus_igh_locus.fasta cow_d_denovo.fasta # search for IGHD genes in the cow IGHD locus
