============
Patch notes
============


v1.0.1
===============

    - Introduced the -u / --genes flag. The output will be limited to the gene IDs/symbols in that file. 
    - Introduced the -z / --reshape flag, which will write the gene-TF matrices in binary format.
    - Updated the TF motif collection (2.2).
    - The conversion of PSCMs to PSEMs will now be done automatically with the average GC-content of the --bed_file as background.
    - One additional peak feature was added to the gene-TF marix, now we report NumPeaks, AvgPeakSize and AvgPeakDistance.
    - The processing order for writing the gene-TF matrices was changed. This will cause a slower calculation for one or very few activity columns, but a faster runtime for more.
    - Fixed an issue that the contact frequency pseudocount was only calculated within the range covered by gene windows. The difference in the resulting pseudocount are negligible.