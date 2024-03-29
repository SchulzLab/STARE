[![Documentation Status](https://readthedocs.org/projects/stare/badge/?version=latest)](https://stare.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5841991.svg)](https://doi.org/10.5281/zenodo.5841991)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/stare-abc/README.html)
[![Platforms](https://anaconda.org/bioconda/stare-abc/badges/platforms.svg)](https://bioconda.github.io/recipes/stare-abc/README.html)


# STARE - Look into TF regulation on gene level leveraging an adapted Activity-By-Contact score

STARE is a framework that combines two tasks, which can also work independently:

 - score enhancer-gene interactions with the Activity-By-Contact (ABC) model
 - derive TF affinities for regions
 
Those two parts are combined to summarise TF affinities on gene level, which allows for a variety of downstream analyses. If this sounds interesting to you, here's a full documentation on [Read the Docs](https://stare.readthedocs.io/en/latest/index.html).

<img src="Figures/STARE_Logo.png" alt="STARE_Logo" width="500"/>

# References

STARE and its application are described in our publication:

Hecker, D., Behjati Ardakani, F., Karollus, A., Gagneur, J., & Schulz, M. H. (2023). The adapted Activity-By-Contact model for enhancer–gene assignment and its application to single-cell data. Bioinformatics, 39(2), btad062. https://doi.org/10.1093/bioinformatics/btad062


The results of our analyses are available via [Zenodo](https://doi.org/10.5281/zenodo.5841991).


The framework is a combination and adaptation which is based on:

- Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). [https://www.nature.com/articles/s41588-019-0538-0](https://www.nature.com/articles/s41588-019-0538-0)
	 - https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
- Combining transcription factor binding affinities with open-chromatin data for accurate gene expression prediction Schmidt et al., Nucleic Acids Research 2016; doi: 10.1093/nar/gkw1061
	 - https://github.com/schulzlab/tepic
- Roider, H. G., A. Kanhere, T. Manke, and M. Vingron. “Predicting Transcription Factor Affinities to DNA from a Biophysical Model.” Bioinformatics 23, no. 2 (January 15, 2007): 134–41. https://doi.org/10.1093/bioinformatics/btl565.


