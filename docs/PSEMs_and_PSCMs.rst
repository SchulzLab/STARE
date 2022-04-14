============
PSCMs and PSEMs
============


There are two possible formats to provide TF motifs to STARE: either Position Specific Count Matrices (PSCMs) in transfac-format or Position Specific Energy Matrices (PSEMs). The more common ones are the PSCMs, which you can get from different public databases. Internally, STARE requires PSEMs to calculate TF affinities with `TRAP <https://doi.org/10.1093/bioinformatics/btl565>`_, and it will automatically convert your PSCMs to PSEMs. The background GC-content to construct the PSEMs will be taken from your provided --bed_file. You can also specify a fixed GC-content with the -y flag. Here's an example for a transfac PSEM from `JASPAR <https://jaspar.genereg.net/>`_::

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

We have multiple motif collections available, which you can find at https://github.com/SchulzLab/STARE/tree/main/PWMs. The newest ones for human and mouse in **2.2/** represent non-redundant motifs collected from:

- `JASPAR 2022 <https://jaspar.genereg.net/>`_: Castro-Mondragon, Jaime A, Rafael Riudavets-Puig, Ieva Rauluseviciute, Roza Berhanu Lemma, Laura Turchi, Romain Blanc-Mathieu, Jeremy Lucas, et al. “JASPAR 2022: The 9th Release of the Open-Access Database of Transcription Factor Binding Profiles.” Nucleic Acids Research 50, no. D1 (January 7, 2022): D165–73. https://doi.org/10.1093/nar/gkab1113.

- `HOCOMOCO v11 <https://hocomoco11.autosome.ru/>`_ : Kulakovskiy, Ivan V, Ilya E Vorontsov, Ivan S Yevshin, Ruslan N Sharipov, Alla D Fedorova, Eugene I Rumynskiy, Yulia A Medvedeva, et al. “HOCOMOCO: Towards a Complete Collection of Transcription Factor Binding Models for Human and Mouse via Large-Scale ChIP-Seq Analysis.” Nucleic Acids Research 46, no. D1 (January 4, 2018): D252–59. https://doi.org/10.1093/nar/gkx1106.

- and from Kheradpour, Pouya, and Manolis Kellis. “Systematic Discovery and Characterization of Regulatory Motifs in ENCODE TF Binding Experiments.” Nucleic Acids Research 42, no. 5 (March 1, 2014): 2976–87. https://doi.org/10.1093/nar/gkt1249.


We have older collections of PSEMs, constructed based on the overall GC-content of the respective organism. In the directory **1.0/** you find the converted PSEMs of vertebrate TFs. In **2.0/** there are two files for each species, one with a collection from different databases *'_all'*, and one where the redundant ones were removed. In **2.1/** there are merged and clustered collections of PSEMs.

If you wish to transform your PSCMs manually with the organism's GC-content as background, here is a selection:

- *Homo sapiens* = 0.41
- *Mus musculus* = 0.42
- *Rattus norvegicus* = 0.42
- *Drosophila melanogaster* = 0.43
- *Caenorhabditis elegans* = 0.36

And the command is the following ::

   ./Code/PSCM_to_PSEM path_to_PSCM GC-content > out_path


For details on TRAP:
- Roider, H. G., A. Kanhere, T. Manke, and M. Vingron. “Predicting Transcription Factor Affinities to DNA from a Biophysical Model.” Bioinformatics 23, no. 2 (January 15, 2007): 134–41. https://doi.org/10.1093/bioinformatics/btl565.


