
Patch notes
============

v1.0.5
===============

 - The STARE_ABCpp binary now creates its own metadata file, listing the used command that called the script along with all the flags.
 - The printouts from STARE_ABCpp when parameters are missing is now more informative. 
 - Gtf-files are now handled more flexibly, catching cases of non-reference chromosomes and missing gene names.
 - Minor fixes for the test cases.
 

v1.0.4
===============

 - From v1.0.3.2: Fixed a bug that caused the Gene-TF-matrices to not include the last batch of genes (max 999).
 - Now allows a **promoter-mode** when leaving the -b flag empty. That means that STARE will build a promoter window of size -w around the 5'TSS of each gene and summarise the affinities in those instead.
 - With >200 activity columns the parallelization of STARE_ABCpp switches from chromosomes to gene batches to allow a better scaling for many columns and to limit the required memory.
 - Catch more cases of combinations of the chr-Prefix in the input files.
 - Now STARE_ABCpp skips chromosomes without genes, e.g. when filtering for genes with -u.
 - Fixed that the bin size -k was required although the contact estimate based on distance was selected.


v1.0.3
===============

 - We found another improvement over the regular ABC score by including information across all annotated TSS of a gene, which we call the generalised ABC (gABC) score. This is the default option now for the ABC part.
 - According to the previous point we added a flag (-i / --tss_mode) which allows to decide whether to use all TSS (all_tss) or just the 5' (5_tss).
 - We added the option to run the ABC-scoring without any chromatin contact data, by setting the -f / --contact_folder to false. The inverse of the distance will then be used as contact estimate.
 - Now properly handling the 1-based gtf-locations when doing intersections with bedtools.


v1.0.2
===============

   - Genes which don't have any regions in their defined window to sum up the TF affinities from will no longer be written to output, only to the discarded_Genes.txt file.
   - Restructured the processing to be more efficient for many activity columns.
   - Added a CMakeList-file to allow compilation with CMake.
   - Window size can now also be given as scientific notation (e.g. 1e05). Before, this could have been an issue if system calls from other programs automatically converted large integers to scientific notation.
   - Fixed a bug that caused unnecessary memory consumption. Further edited some data structures to reduce memory usage.
   - Fixed a system-dependent bug when the limit of open streams was reached. Out-streams are now always closed after writing.
   - Fixed a bug related to the thread-safety of streams.


v1.0.1
===============

    - Introduced the -u / --genes flag. The output will be limited to the gene IDs/symbols in that file. 
    - Introduced the -z / --reshape flag, which will write the gene-TF matrices in binary format.
    - Updated the TF motif collection (2.2).
    - The conversion of PSCMs to PSEMs will now be done automatically with the average GC-content of the --bed_file as background.
    - One additional peak feature was added to the gene-TF marix, now we report NumPeaks, AvgPeakSize and AvgPeakDistance.
    - The processing order for writing the gene-TF matrices was changed. This will cause a slower calculation for one or very few activity columns, but a faster runtime for more.
    - Fixed an issue that the contact frequency pseudocount was only calculated within the range covered by gene windows. The difference in the resulting pseudocount are negligible.