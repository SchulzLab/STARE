#!/bin/bash -e

# Part of STARE: https://github.com/SchulzLab/STARE
# Small bash script to ease the normalization of .hic-files with Juicebox (https://github.com/aidenlab/juicer/wiki/Data-Extraction)

help="\n
Usage: ./Juicebox_KR_normalization.sh [-h hic-file to normalize]\n[-j path to the jar-file]\n
[-d folder to write the normalized files to]\n
optional:\n
[-c the chromosomes to normalize, e.g. 1-22+XY, default is 1-22, also possible to query individual chromosomes or multiple ones comma separated (1,5,7)]\n
[-b bin size, meaning resolution of the matrix (default 5000)]\n"

# ------------------------------------------------------------------------------------------------------
# FETCHING INPUT
# ------------------------------------------------------------------------------------------------------
hic_file=""
out_folder=""
jar_file=""
chromosomes="1-22"
bin_size=5000

# Parsing command line.
while getopts "h:d:j:c:b:" o;
do
case $o in
	h)	hic_file=$OPTARG;;
	d)	out_folder=$OPTARG;;
	j)	jar_file=$OPTARG;;
	c)	chromosomes=$OPTARG;;
  b)  bin_size=$OPTARG;;
esac
done

if [ -z "$hic_file" ] ;
then
	echo Hi-C file must be specified with -h
	exit 1;
fi

if [ -z "$out_folder" ] ;
then
	echo An output folder must be specified with -d
	exit 1;
fi

if [ -z "$jar_file" ] ;
then
	echo The jar file is needed to call Juicebox
	exit 1;
fi

# Create the out_dir if not existent and call Juicebox for each of the specified chromosomes.
if [ ! -d "$out_folder" ]; then
  mkdir "$out_folder"
fi

if [[ $chromosomes == *"-"* ]]; then
  chr_range=$(echo ${chromosomes} | tr -d 'XY')
  IFS='-' read -ra ADDR <<< "${chr_range}"
  for i in $(seq "${ADDR[0]}" "${ADDR[1]}"); do
    java -jar ${jar_file} dump observed KR "${hic_file}" chr${i} chr${i} BP ${bin_size} ${out_folder}/chr${i}_KR_Contacts.txt;
    gzip ${out_folder}/chr${i}_KR_Contacts.txt;
  done
  # Manually check if any gonosome was selected in addition to the range.
  if [[ $chromosomes == *"X"* ]]; then
    java -jar ${jar_file} dump observed KR ${hic_file} chrX chrX BP ${bin_size} ${out_folder}/chrX_KR_Contacts.txt;
    gzip ${out_folder}/chrX_KR_Contacts.txt;
  fi
  if [[ $chromosomes == *"Y"* ]]; then
    java -jar ${jar_file} dump observed KR ${hic_file} chrY chrY BP ${bin_size} ${out_folder}/chrY_KR_Contacts.txt;
    gzip ${out_folder}/chrY_KR_Contacts.txt;
  fi
elif [[ $chromosomes == *","* ]]; then
    IFS=',' read -ra csv <<< "$chromosomes"
    for i in "${csv[@]}"; do
    java -jar ${jar_file} dump observed KR ${hic_file} chr${i} chr${i} BP ${bin_size} ${out_folder}/chr${i}_KR_Contacts.txt;
    gzip ${out_folder}/chr${i}_KR_Contacts.txt;
    done
else  # If only one individual chromosome was queried.
    java -jar ${jar_file} dump observed KR "${hic_file}" chr${chromosomes} chr${chromosomes} BP ${bin_size} ${out_folder}/chr${chromosomes}_KR_Contacts.txt;
    gzip ${out_folder}/chr${chromosomes}_KR_Contacts.txt;
fi

