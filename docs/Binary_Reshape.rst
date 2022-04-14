============
Binary output
============

You have the option to receive a binary output file from STARE with the -z flag, which replaces the default human readable files. This should mostly be interesting for you if you run STARE together with `GAZE <https://github.com/schulzlab/gaze>`_. If you would like to have the binary version independent of GAZE, for example for faster processing, here is an explanation of the format.

*ExampleOutput_Reshape_Binary.bin.gz* and *ExampleOutput_Reshape_Meta.txt*
===============
These two files together contain your data: *ExampleOutput_Reshape_Binary.bin.gz* the floats and *ExampleOutput_Reshape_Meta.txt* the accompanying metadata. The floats are represented in IEEE-754 format with 32-bit precision, meaning that each block of 32 bits is one float. The byte order is dependent on the system where the file was written to, see the section below. If the binary format was selected, the output of all activity columns are written into one file. The order is TFs-column-gene, as depicted below::

    TF1-Column1-Gene1 | TF2-Column1-Gene1 ... TFn-Column1-Gene1 | TF1-Column2-Gene1...

That means that if you analysed 100 TFs, the first 100 entries in the binary files will be the scores for column1 for gene1, followed by 100 entries for column2 for gene1. The entries for Gene2 start at *index number of TFs * number of columns*. The order of TFs, activity columns and genes is listed in the *ExampleOutput_Reshape_Meta.txt* file. The first row in the meta file however, is the number of floats stored in the binary file. That means that the byte size of the uncompressed binary file should be the number of floats times 4.

Byte-order / Endianness
===============
As uncomplicated interoperability across systems would be boring, there are different orders to store bytes in a computer memory. This byte order, the so called endianness, decides how the floats in our binary file are stored. Forcing the byte-order upon file creation turned out not to be straightforward, so we provide two additional files which should help in identifying the endianness of the system on which the binary file was written. *_Endianness_check.txt* names the endianness in the first row, which should either be LITTLE_ENDIAN or BIG_ENDIAN. The next lines in this file are float numbers. Those float numbers are also stored in *_Endianness_check.bin*, but in their binary representation (also 32-bit IEEE-754). If you like, you can use those to check whether the read-in works properly and whether you get the same numbers from the bin- as from the txt-file. We chose floats which should be unambiguous to any rounding rules.

