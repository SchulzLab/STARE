============
Separate ABC-scoring
============

If you're only interested in the ABC-scored interactions, you can also call that part independently. The flags are the same as for whole STARE, although the long options are not available. Here is an example on how to run it::

   ./Code/STARE_ABCpp -b <path_to_bed_file> -n <activity_column(s) -a <gtf_annotation> -o <output_path> -w <window_size> -f <contact_data_dir> -k <bin_size> -t <score_cut_off>


For convenience, here are the flags that are specific to the ABC-scoring alone:

**required**

.. csv-table:: 
   :header: "Flag", "Description"
   :widths: 8, 50

   -b, Bed-file containing your candidate regions. Headers are allowed if they start with #.
   -a, Gene annotation file in gtf-format.
   -o, Output-prefix with which all file names will start. Unlike whole STARE no separate folder will be created.

**optional**

.. csv-table:: 
   :header: "Flag", "Description"
   :widths: 8, 50

   -n, Column(s) in the --bed_file representing the activity of the region. You will get one set of output files for each column. Start counting at 1. Allowed formats are individual columns; column ranges; columns separated by comma as well as a start column with all consecutive columns.  
   -c, Number of cores to provide for parallel computing. Note that the processing is also heavy on memory.
   -x, Bed-file with regions to exclude. All regions in the --bed_file with ≥ 1 bp overlap will be discarded from all further analyses. Make sure that your -b file has the same naming scheme e.g. both *chr1* or both *1* in the first column. Otherwise an intersection is not possible.
   -w, Window size centred at the 5' TSS in which regions from the --bed_file will be considered for a gene (Default 50KB for non-ABC-mode and 5MB for ABC-mode). E.g. 5MB means ±2.5MB around the TSS.
   -f, Path to directory containing normalized chromatin contact files in coordinate format (bin|bin|contact) one gzipped file for each chromosome.
   -k,  Resolution of the chromatin contact data. E.g. 5000 for a 5kb resolution.
   -t,  Cut-off for the ABC-score. Only interactions surpassing it are written to the output (Default 0.02). Set to 0 if you would like to get all scored interactions.
   -q,  Whether to use the use the adapted ABC-scoring or the 'original' one (Default True).
   -m,  Size of the window around your candidate regions in which genes are considered for the adapted activity adjustment (Default 5MB; will be minimally set to -w).
   -d,  Whether to use a pseudocount for the contact frequency in the ABC-score (Default True).
   -h, Print the flag options.
