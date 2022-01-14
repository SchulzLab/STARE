# STARE - Look into TF regulation on gene level leveraging an adapted Activity-By-Contact score

STARE is a framework that combines two tasks, which also work independently from each other:

 - score enhancer-gene interactions with the Activity-By-Contact (ABC) model
 - derive TF affinities for regions
 
Those two parts are combined to summarise TF affinities on gene level, which allows for a variety of further downstream analyses. If this sounds interesting to you, here's how you get started.

# Installation
*We are working on bringing STARE to Bioconda. But until then, you need to do the installation manually.*  
* [bedtools](https://github.com/arq5x/bedtools2); Please make sure to add the bedtools installation to your path
* g++ to compile the C++ scripts 
* openmp for parallel computing; unfortunately, the installation is system-dependent

Once you got all of the above installed, you can find two scripts in **/Code**: *compile_STARE_macOS.sh* and *compile_STARE_Linux.sh* which should compile the C++ scripts for your respective system if you run the following command from within the **/Code** directory:

    ./compile_STARE_macOS.sh
    <or>
    ./compile_STARE_Linux.sh

# Get started
The following schema should give you an overview on how you can tune STARE's settings to run it properly for your data. 

![STARE_Flow](/Figures/STARE_FlowBig.png)

If you want to test your installation and try out some example runs, look into the **/Test** folder. From within the folder you can run:

    ./run runTestCases.sh

It executes a range of very small examples with varying flags and should tell you if something didn't go as planned.

# Output
The output unsurprisingly depends on your input and the options you chose. However, the overall structure is similar, and you might not produce all files that are listed here. Let's pretend you set the -o flag to *pancake* so that we can have the full paths as example.

 You will always get a metadata file *Pancake_metadata.amd.tsv*, which lists the flags you set.

**ABC output:**

- *Pancake_ABCpp_scoredInteractions.txt.gz* : Lists all the enhancer-gene interactions that were scored higher than your set cut-off -t. 
 - *Pancake_GeneInfo.txt.gz*: Summarises a variety of features for each gene. Here are details on the not so intuitive ones:
	 - *Avg_EnhancerActivity*: average of the activity of all enhancers that were mapped to that gene
	 - *Avg_EnhancerContact* average contact of all enhancers that were mapped to that gene, that includes the pseudocount
	 - *Avg_EnhancerDistance*: average distance of the enhancers to the TSS of the gene
	 - *Failure*: in case that the gene could not been processed or when no enhancer was located within the gene window, you will find a respective note here, otherwise it is '-'

**Gene-TF-matrices:**

 - *Pancake_TF_Gene_Affinities.txt.gz*: Matrix of TF affinities summarised per gene, with the genes as rows and TFs as columns. It has two additional columns with the average peak size and average peak distance of the regions that were considered for the gene.
 - *Pancake_discarded_Genes.txt*: Lists all genes where no TF affinities could be calculated, with a note indicating why.

If you selected multiple activity columns (-n), you will also receive output files for each.


# References
We will soon provide a preprint paper about STARE and its implementation, including application examples. Up to this point, note that STARE is a combination and adaptation which is based on:

 - Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). [https://www.nature.com/articles/s41588-019-0538-0](https://www.nature.com/articles/s41588-019-0538-0)
	 - https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
 - Combining transcription factor binding affinities with open-chromatin data for accurate gene expression prediction Schmidt et al., Nucleic Acids Research 2016; doi: 10.1093/nar/gkw1061
	 - https://github.com/schulzlab/tepic
