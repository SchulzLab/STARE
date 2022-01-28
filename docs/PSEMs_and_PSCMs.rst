============
PSEMs and PSCMs
============

The TF affinities are calculated with `TRAP <https://doi.org/10.1093/bioinformatics/btl565>`_. TRAP is motif based and requires Position Specific Energy Matrices (PSEMs). They represent the mismatch energy of a given motif. The authors of TRAP provide a tool which converts Position Specific Count Matrices (PSCMs) to PSEMs. We provide a collection of PSEMs for different species, which you can find at https://github.com/SchulzLab/STARE/tree/main/PWMs. In the directory **1.0/** you find the converted PSEMs of vertebrate TFs. In **2.0** there are two files for each species, one with a collection from different databases *'_all'*, and one where the redundant ones were removed. In **2.1/** there are merged and clustered collections of PSEMs. The original PSCMs were collected from:

- `JASPAR <https://jaspar.genereg.net/>`_,: Sandelin, A. “JASPAR: An Open-Access Database for Eukaryotic Transcription Factor Binding Profiles.” Nucleic Acids Research 32, no. 90001 (January 1, 2004): 91D – 94. https://doi.org/10.1093/nar/gkh012.

- `HOCOMOCO <https://hocomoco11.autosome.ru/>`_ : Kulakovskiy, Ivan V., Yulia A. Medvedeva, Ulf Schaefer, Artem S. Kasianov, Ilya E. Vorontsov, Vladimir B. Bajic, and Vsevolod J. Makeev. “HOCOMOCO: A Comprehensive Collection of Human Transcription Factor Binding Sites Models.” Nucleic Acids Research 41, no. D1 (January 1, 2013): D195–202. https://doi.org/10.1093/nar/gks1089.

- and from Kheradpour, Pouya, and Manolis Kellis. “Systematic Discovery and Characterization of Regulatory Motifs in ENCODE TF Binding Experiments.” Nucleic Acids Research 42, no. 5 (March 1, 2014): 2976–87. https://doi.org/10.1093/nar/gkt1249.

If your TF of interest are not represented in any of the provided collections, you can hand STARE your PSCMs and it will automatically convert them. The PSCMs need to be in transfac-format, an example from `JASPAR <https://jaspar.genereg.net/>`_::

	AC MA0006.1
	XX
	ID Ahr::Arnt
	XX
	DE MA0006.1 Ahr::Arnt ; From JASPAR
	PO	A	C	G	T
	01	3.0	8.0	2.0	11.0
	02	0.0	0.0	23.0	1.0
	03	0.0	23.0	0.0	1.0
	04	0.0	0.0	23.0	1.0
	05	0.0	0.0	0.0	24.0
	06	0.0	0.0	24.0	0.0
	XX
	CC tax_group:vertebrates
	CC tf_family:PAS domain factors; PAS domain factors
	CC tf_class:Basic helix-loop-helix factors (bHLH); Basic helix-loop-helix factors (bHLH)
	CC pubmed_ids:7592839
	CC uniprot_ids:P30561; P53762
	CC data_type:SELEX
	XX
	//

The conversion also requires the average GC-content of the organism. Here a selection:

- *Homo sapiens* = 0.41
- *Mus musculus* = 0.42
- *Rattus norvegicus* = 0.42
- *Drosophila melanogaster* = 0.43
- *Caenorhabditis elegans* = 0.36

You can also convert the PSCMs in advance ::

   ./Code/PSCM_to_PSEM path_to_PSCM GC-content > out_path


For details on TRAP:
- Roider, H. G., A. Kanhere, T. Manke, and M. Vingron. “Predicting Transcription Factor Affinities to DNA from a Biophysical Model.” Bioinformatics 23, no. 2 (January 15, 2007): 134–41. https://doi.org/10.1093/bioinformatics/btl565.


