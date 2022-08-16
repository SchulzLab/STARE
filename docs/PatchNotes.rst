
Patch notes
============

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