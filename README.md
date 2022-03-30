[![Documentation Status](https://readthedocs.org/projects/stare/badge/?version=latest)](https://stare.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5841991.svg)](https://doi.org/10.5281/zenodo.5841991)


# STARE - Look into TF regulation on gene level leveraging an adapted Activity-By-Contact score

STARE is a framework that combines two tasks, which can also work independently:

 - score enhancer-gene interactions with the Activity-By-Contact (ABC) model
 - derive TF affinities for regions
 
Those two parts are combined to summarise TF affinities on gene level, which allows for a variety of downstream analyses. If this sounds interesting to you, here's a full documentation on [Read the Docs](https://stare.readthedocs.io/en/latest/index.html).

<img src="Figures/STARE_Logo.png" alt="STARE_Logo" width="500"/>

# References

STARE and its application are described in our preprint on bioRxiv:

Hecker, Dennis, Fatemeh Behjati Ardakani, and Marcel H. Schulz. “The Adapted Activity-By-Contact Model for Enhancer-Gene Assignment and Its Application to Single-Cell Data” Preprint. BioRxiv, January 28, 2022. https://doi.org/10.1101/2022.01.28.478202.


The results of our analyses are available via [Zenodo](https://doi.org/10.5281/zenodo.5841991).


The framework is a combination and adaptation which is based on:

- Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). [https://www.nature.com/articles/s41588-019-0538-0](https://www.nature.com/articles/s41588-019-0538-0)
	 - https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
- Combining transcription factor binding affinities with open-chromatin data for accurate gene expression prediction Schmidt et al., Nucleic Acids Research 2016; doi: 10.1093/nar/gkw1061
	 - https://github.com/schulzlab/tepic
- Roider, H. G., A. Kanhere, T. Manke, and M. Vingron. “Predicting Transcription Factor Affinities to DNA from a Biophysical Model.” Bioinformatics 23, no. 2 (January 15, 2007): 134–41. https://doi.org/10.1093/bioinformatics/btl565.


