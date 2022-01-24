

# STARE - Look into TF regulation on gene level leveraging an adapted Activity-By-Contact score

STARE is a framework that combines two tasks, which can also work independently:

 - score enhancer-gene interactions with the Activity-By-Contact (ABC) model
 - derive TF affinities for regions
 
Those two parts are combined to summarise TF affinities on gene level, which allows for a variety of further downstream analyses. If this sounds interesting to you, here's how you get started.

# Installation
*We are working on bringing STARE to Bioconda. But until then, you need to do the installation manually.*  

* [bedtools](https://github.com/arq5x/bedtools2); Please make sure to add the bedtools installation to your path
* g++ to compile the C++ scripts 
* openmp for parallel computing; unfortunately, the installation is system-dependent
* [Boost C++ Library](https://www.boost.org/)

Once you got all of the above installed, you can clone the GitHub repository

    git clone https://github.com/SchulzLab/STARE.git
or download the source code and unzip it.  You will find two scripts in **/Code**: *compile_STARE_macOS.sh* and *compile_STARE_Linux.sh* which should compile the C++ scripts for your platform if you run them. The paths are relative, you can call them from any directory:

    ./Code/compile_STARE_macOS.sh
or

    ./Code/compile_STARE_Linux.sh

# Get started
The following schema should give you an overview on how you can tune STARE's settings to run it properly for your data. 

![STARE_Flow](/Figures/STARE_FlowBig.png)

If you want to test your installation and try out some example runs, we have the *Code/runTestCases.sh* script for you. It serves the following purposes: 
- It gives examples on how to run STARE and which flags to use. To get inspiration have a look at the individual test commands.
- It also compares the output of its test runs against pre-computed results to make sure that the installation worked correctly (**/Test/Test_Controls/**). It will tell you with a subtle ERROR massage if something went wrong.

The paths are relative, you should be able to call the test suite from anywhere:

    ./Code/runTestCases.sh


# Output
The output depends on your input and the options you chose. Although you might not produce all files that are listed, the overall structure is similar. Let's pretend you set the -o flag to *pancake* so that we can have the full paths as example.

 You will always get a **metadata** file *Pancake_metadata.amd.tsv*, which lists the flags you set.

**ABC output:**

- *Pancake_ABCpp_scoredInteractions.txt.gz* : Lists all the enhancer-gene interactions that were scored higher than your set cut-off -t. 
 - *Pancake_GeneInfo.txt.gz*: Summarises a variety of features for each gene. Here are details on the not so intuitive ones:
	 - *Avg_EnhancerActivity*: average of the activity of all enhancers that were mapped to that gene
	 - *Avg_EnhancerContact*: average contact of all enhancers that were mapped to that gene, that includes the pseudocount
	 - *Avg_EnhancerDistance*: average distance of the enhancers to the TSS of the gene
	 - *Failure*: in case that the gene could not be processed or when no enhancer was located within the gene window, you will find a respective note here, otherwise it is '-'

**Gene-TF-matrices:**

 - *Pancake_TF_Gene_Affinities.txt.gz*: Matrix of TF affinities summarised per gene, with the genes as rows and TFs as columns. It has two additional columns with the average peak size and average peak distance of the regions that were considered for the gene.
 - *Pancake_discarded_Genes.txt*: Lists all genes where no TF affinities could be calculated, with a note indicating why.

If you selected multiple activity columns (-n), you will also receive output files for each. Start counting at 1. You can indicate the columns in multiple ways. Here some examples:

 - **-n 4**  :  use the activity in column 4
 - **-n 4-6** : use column 4, 5 and 6
 - **-n 4+** : run all consecutive columns, from column 4 to the last one




# References
We will soon provide a preprint paper about STARE and its implementation, including application examples. Up to this point, note that STARE is a combination and adaptation which is based on:

 - Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). [https://www.nature.com/articles/s41588-019-0538-0](https://www.nature.com/articles/s41588-019-0538-0)
	 - https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
 - Combining transcription factor binding affinities with open-chromatin data for accurate gene expression prediction Schmidt et al., Nucleic Acids Research 2016; doi: 10.1093/nar/gkw1061
	 - https://github.com/schulzlab/tepic

The results of our analyses are available via Zenodo.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5841992.svg)](https://doi.org/10.5281/zenodo.5841991)
