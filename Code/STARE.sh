#!/bin/bash -e
set -e  # To abort the whole script if one function returns an error.

# See https://github.com/SchulzLab/STARE for more information and usage.
# Adapted from TEPIC: https://github.com/SchulzLab/TEPIC
help="STARE version 0.1\n
Usage: ./STARE.sh
[-b bed file containing open chromatin regions]\n
[-g input fasta file in RefSeq format]\n
[-j flag indicating that the reference genome contains a chr prefix]\n
[-p file with PSEMS of TFs]\n
[-a gene annotation file, required to generate the gene view]\n
[-o prefix_path of output files]\n
Optional parameters:\n
[-n column in the -b file containing the average per base signal within a peak, start counting at 1]\n
[-c number of cores to use (default 1)]
[-x bed-file with regions to exclude (e.g. blacklisted regions)]\n
[-w window size around TSS for mapping regions to genes (default 50KB; 5MB for ABC-mode)]\n
[-e indicating whether exponential distance decay should be used (default TRUE, but not used in ABC-mode)]\n
[-f folder with normalized Hi-C contact files for each chromosome in coordinate format]\n
[-k bin-size of the Hi-C files]\n
[-t cut-off for the ABC-score (default 0.02), set to 0 to get all scored interactions]\n
[-q whether to use the adapted ABC-score (default True)]\n
[-m window size around enhancers for the -q adjustment (default 5MB, minimally set to -w)]\n
[-d whether to use pseudocount for the contact frequency in the ABC-model (default True)]\n
[-r ABC-scoring file, if already calculated once for this input to avoid redundant calculation]\n
\n"

# ------------------------------------------------------------------------------------------------------
# FETCHING INPUT
# ------------------------------------------------------------------------------------------------------
genome=""
regions=""
prefixP=""
cores=1
pwms=""
column="0"
annotation=""
window=""
decay="TRUE"
chrPrefix="FALSE"
hic_contactfolder=""
hic_binsize=""
abc_cutoff=0.02
existing_abc="0"
exclude_regions=""
pseudocount="TRUE"
adjustedABC="TRUE"
enhancer_window=5000000

# Parsing command line.
while getopts "g:b:o:c:p:n:a:w:e:j:f:k:t:r:x:d:q:m:" o;
do
case $o in
	g)	genome=$OPTARG;;
	b)	regions=$OPTARG;;
	o)	prefixP=$OPTARG;;
	c)	cores=$OPTARG;;
	p)	pwms=$OPTARG;;
	n)	column=$OPTARG;;
	a)	annotation=$OPTARG;;
	w)	window=$OPTARG;;
	e)	decay=$OPTARG;;
	j)	chrPrefix=$OPTARG;;
  f)  hic_contactfolder=$OPTARG;;
  k)  hic_binsize=$OPTARG;;
  t)  abc_cutoff=$OPTARG;;
  r)  existing_abc=$OPTARG;;
  x)  exclude_regions=$OPTARG;;
  d)  pseudocount=$OPTARG;;
  q)  adjustedABC=$OPTARG;;
  m)  enhancer_window=$OPTARG;;
  *)  exit 1;;
esac
done

if [ $OPTIND -eq 1 ] ;
then
    echo -e "$help"
    exit 1;
fi

if [ -z "$genome" ] ;
then
	echo Reference genome must be specified using the -g parameter
	exit 1;
fi

if [ -z "$annotation" ] ;
then
	echo Gene annotation must be specified using the -g parameter
	exit 1;
fi

if [ -z "$regions" ] ;
then
	echo Open chromatin regions must be specified using the -b parameter
	exit 1;
fi

if [ -z "$prefixP" ] ;
then
	echo Prefix of output files must be specified using the -o parameter
	exit 1;
fi

if [ -z "$pwms" ] ;
then
	echo PWMs must be specified using the -p parameter
	exit 1;
fi

if [ -n "$hic_contactfolder" ] || [ -n "$hic_binsize" ] && [[ "$existing_abc" == "0" ]];  # With an existing ABC-file, the other flags are ignored.
then
  if [ -z "$hic_contactfolder" ] || [ -z "$hic_binsize" ] || [ -z "$column" ];
  then
    echo "For the ABC-score calculation the column with the peak signal (-n), the path to the normalized contact files (-f) as well as the the bin size (-k) are needed."
    exit 1;
  fi
fi

if [ -z "$window" ] ;
then
  window=50000
  if [ -n "$hic_contactfolder" ] || [ -n "$hic_binsize" ] ;  # Change default if the ABC-mode is used.
  then
    window=5000000
  fi
fi

if [ "${chrPrefix}" == "TRUE" ] || [ "${chrPrefix}" == "True" ] || [ "${chrPrefix}" == "true" ] || [ "${chrPrefix}" == "T" ] || [ "${chrPrefix}" == "1" ] ;
then
  chrPrefix="TRUE";
fi

if [ "${decay}" == "FALSE" ] || [ "${decay}" == "False" ] || [ "${decay}" == "false" ] || [ "${decay}" == "F" ] || [ "${decay}" == "0" ] ;
then
  decay="FALSE";
fi

if [ "${pseudocount}" == "FALSE" ] || [ "${pseudocount}" == "False" ] || [ "${pseudocount}" == "false" ] || [ "${pseudocount}" == "F" ] || [ "${pseudocount}" == "0" ] ;
then
  pseudocount="FALSE";
fi

if [ "${adjustedABC}" == "FALSE" ] || [ "${adjustedABC}" == "False" ] || [ "${adjustedABC}" == "false" ] || [ "${adjustedABC}" == "F" ] || [ "${adjustedABC}" == "0" ] ;
then
  adjustedABC="FALSE";
fi

# ------------------------------------------------------------------------------------------------------
# WRITE METADATA FILE
# ------------------------------------------------------------------------------------------------------
d=$(date +%D)
d=`echo $d | sed 's/\//\_/g'`
t=$(date +%H:%M:%S:%N | sed 's/:/_/g')

if [ -d "$prefixP" ]; then
  echo "$prefixP" " Folder already exists, remove it or change the function call prefix in the -o flag"
  exit 1;
fi

mkdir "${prefixP}"
base_prefix=$(basename "${prefixP}")
prefix_path=$prefixP"/"$base_prefix
working_dir=$(cd "$(dirname "$0")" && pwd -P)

filteredRegions=$prefix_path"_candiate_binding_regions"
# Generating name of the fasta file containing the overlapping regions.
openRegionSequences=${prefix_path}.OpenChromatin.fasta
metadatafile=${prefix_path}_metadata.amd.tsv
# Create metadata file.
touch "$metadatafile"
echo "[Description]" >> "$metadatafile"
echo "process	STARE0.1" >> "$metadatafile"
echo -e "run_by_user\t""$USER" >> "$metadatafile"
echo -e "date\t""$d" >> "$metadatafile"
echo -e "time\t""$t" >> "$metadatafile"
echo -e "analysis_id\t""$prefix_path" >> "$metadatafile"
echo "" >> "$metadatafile"
echo "[Inputs]" >> "$metadatafile"
echo -e "region_file\t""$regions" >> "$metadatafile"
if [ -n "$column" ] ;
then
	echo -e "signal_column\t""$column" >> "$metadatafile"
fi
if [ -n "$exclude_regions" ] ;
then
	echo -e "excluded regions\t""$exclude_regions" >> "$metadatafile"
fi
echo "" >> "$metadatafile"
echo "[References]" >> "$metadatafile"
echo -e "genome_reference\t""$genome" >> "$metadatafile"
echo -e "pwms\t""$pwms" >> "$metadatafile"
echo -e "genome_annotation\t""$annotation">> "$metadatafile"

echo "" >> "$metadatafile"
echo "[Output path]" >> "$metadatafile"
echo "$prefixP" >> "$metadatafile"

echo "" >> "$metadatafile"
echo "[Parameters]" >> "$metadatafile"
echo -e "SampleID\t""$prefixP" >> "$metadatafile"
echo -e "cores\t""$cores" >> "$metadatafile"
echo -e "chr prefix\t"$chrPrefix >> "$metadatafile"
echo -e "window\t"$window >> "$metadatafile"
if [ -z "$hic_contactfolder" ] && [[ "$existing_abc" == "0" ]];
then
  echo -e "decay\t""$decay" >> "$metadatafile"
fi
if [ -n "$hic_contactfolder" ] && [[ "$existing_abc" == "0" ]];
then
  echo -e "path with hi-c contact files\t""$hic_contactfolder" >> "$metadatafile"
  echo -e "bin size of hi-c contacts\t""$hic_binsize" >> "$metadatafile"
  echo -e "ABC-score cut-off\t""$abc_cutoff" >> "$metadatafile"
  echo -e "Use pseudocount for contact frequency\t""$pseudocount" >> "$metadatafile"
  echo -e "Use adjustedABC version\t""$adjustedABC" >> "$metadatafile"
  echo -e "Window size for the adjustedABC\t""$enhancer_window" >> "$metadatafile"
fi
if [[ "$existing_abc" != "0" ]];
then
  echo -e "existing ABC-score file that was used\t""$existing_abc" >> "$metadatafile"
fi
echo "" >> "$metadatafile"
echo "[Metrics]" >> "$metadatafile"
numReg=`grep -c . "$regions"`
echo -e "Number of provided regions\t""$numReg" >> "$metadatafile"
numMat=`grep ">" "$pwms" | wc -l`
echo -e "Number of considered pwms\t""$numMat" >> "$metadatafile"

# ------------------------------------------------------------------------------------------------------
# REGION PROCESSING
# ------------------------------------------------------------------------------------------------------
echo "Preprocessing region file"

sed 's/chr//g' "$regions" >  "${filteredRegions}"_Filtered_Regions.bed
sort -s -k1,1 -k2,2 -k3,3 "${filteredRegions}"_Filtered_Regions.bed | uniq > "${filteredRegions}"_sorted.bed
rm "${filteredRegions}"_Filtered_Regions.bed

# Remove regions that overlap with regions in the $exclude_regions file with the intersect -v flag.
if [ -n "$exclude_regions" ] ;
then
  sed 's/chr//g' "$exclude_regions" >  "${exclude_regions}"_noPrefix.bed
  bedtools intersect -a "${filteredRegions}"_sorted.bed -b "${exclude_regions}"_noPrefix.bed -v -header > "${filteredRegions}"_filtered.bed
  rm "${filteredRegions}"_sorted.bed
  rm "${exclude_regions}"_noPrefix.bed
  mv "${filteredRegions}"_filtered.bed "${filteredRegions}"_sorted.bed  # Rename to the original file again.
fi

if [ "${chrPrefix}" == "TRUE" ];
then
	awk '{print "chr"$1"\t"$2"\t"$3}' "${filteredRegions}"_sorted.bed > "${prefix_path}"_tempFasta.bed
	sed -e '/^chr#/d' "${prefix_path}"_tempFasta.bed > "${prefix_path}"_tempFastaNoHeader.bed
	rm "${prefix_path}"_tempFasta.bed
	getFastaRegion=${prefix_path}_tempFastaNoHeader.bed
else
	getFastaRegion=${filteredRegions}_sorted.bed
fi

#echo "Running bedtools"
# Run bedtools to get a fasta file containing the sequence data for predicted open chromatin regions contained in the bedfile.
bedtools getfasta -fi "$genome" -bed "${getFastaRegion}" -fo "$openRegionSequences"
if [ "${chrPrefix}" == "TRUE" ];
then
	rm "${getFastaRegion}"
fi

# Replace invalid characters in the fasta file with 'N', but skip the id rows (">...").
"${working_dir}"/ReplaceInvalidChars -i "$openRegionSequences" -o "${prefix_path}"_FilteredSequences.fa -d "${prefix_path}"_maxRow.txt
rm "$openRegionSequences"

# ------------------------------------------------------------------------------------------------------
# TF-AFFINITY WITH TRAP
# ------------------------------------------------------------------------------------------------------
startt=`date +%s`
# Use TRAP to compute transcription factor affinities to the above extracted sequences.
affinity=${prefix_path}_Affinity.txt

echo "Starting TRAP"
"${working_dir}"/TRAPmulti "$pwms" "${prefix_path}"_FilteredSequences.fa "${prefix_path}"_maxRow.txt "$cores" > "${affinity}"
rm "${prefix_path}"_FilteredSequences.fa
rm "${prefix_path}"_maxRow.txt

endt=`date +%s`
echo $((endt-startt))"s TRAP"

# ------------------------------------------------------------------------------------------------------
# RUN ABC-SCORING if flags are present
# ------------------------------------------------------------------------------------------------------

if [ -n "$hic_contactfolder" ] && [[ "$existing_abc" == "0" ]];  # This point isn't reached if the bin size isn't set as well.
then
  echo "ABC-scoring region-gene interactions"
  mkdir "${prefixP}""/ABC_output"
  abc_prefix_path=${prefixP}"/ABC_output/"${base_prefix}
  "${working_dir}"/STARE_ABCpp -b "${filteredRegions}"_sorted.bed -n "${column}" -a "${annotation}" -gw "${window}" -cf "${hic_contactfolder}" -bin "${hic_binsize}" -t "${abc_cutoff}" -d "${abc_prefix_path}" -p "${pseudocount}" -q "${adjustedABC}" -m "${enhancer_window}" -c "${cores}"
  existing_abc=${abc_prefix_path}"_ABCpp_scoredInteractions.txt.gz"
fi

# ------------------------------------------------------------------------------------------------------
# GET TF-GENE AFFINITIES
# ------------------------------------------------------------------------------------------------------
# The gene view is generated.
startg=`date +%s`
echo "Generating TF-Gene scores"
mkdir "${prefixP}""/Gene_TF_matrices"
"${working_dir}"/TF_Gene_Scorer -a "${annotation}" -b "${filteredRegions}"_sorted.bed -n "${column}" -i "${affinity}" -o "${prefixP}"/Gene_TF_matrices/"${base_prefix}" -p "${pwms}" -w ${window} -e "${decay}" -c "${cores}" -abc "${existing_abc}"

endg=`date +%s`
echo $((endg-startg))"s TF-Gene Scores"


# Clean-Up
rm "${affinity}"
rm "${filteredRegions}"_sorted.bed

echo "Congratulations it worked!"
