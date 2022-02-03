============
Troubleshooting
============

We made STARE rather verbose, so that it tells you what it is currently doing. Also we try to catch the more obvious issues, like unavailable files, and print an appropriate error message. Nonetheless, we can't cover all possibilities, especially since there are so many different data types coming together. We thought of some things that might cause issues and that might be worth checking, if you encounter an error. The list is not exhaustive, it is more driven by personal experience.

Index of activity columns
******************
Make sure that all the columns you gave with the -n flags actually exist and that you start counting at 1. Otherwise it might cause an index error, C++ will call a segmentation fault.

Header in bed-files
******************
Header are allowed in bed-files, and also encouraged to label the activity columns, but they must be preceeded by a number sign '#'. Otherwise, bedtools will fail to work.

The *chr* prefix
******************
We rely on `bedtools <https://bedtools.readthedocs.io/en/latest/>`_ for intersection of regions. Bedtools requires that the naming scheme for the first column is identical, meaning that either both files have a chr-prefix, or both don't have it. A mixture will cause an error. STARE will check the files in advance and should take care of this for you. But if there is an unexpected format or the scaffolds don't start with "chr", this will likely fail. This affects the region file (-b), gene annotation (-a), sequence file (-g) and regions to exclude (-x). It might be worth ensuring consistency for the chromosome name in each of the files, if you suspect that to be the issue.

Is it compressed?
******************
Almost all files are gzipped by STARE. If you decompress a file and remove the compressed version, STARE won't find the file it needs.

Still not working?
******************
If you still have trouble getting STARE running, open an `issue on GitHub <https://github.com/SchulzLab/STARE/issues>`_ and we can try to solve it together. Also don't hesitate to open an issue if you have suggestions, feedback or ideas for STARE.
