#!/bin/bash -e
set -e  # To abort the whole script if one function returns an error.

# See https://github.com/SchulzLab/STARE for more information and usage.
# Adapted from TEPIC: https://github.com/SchulzLab/TEPIC
version_num="1.0.4"
help="STARE version ""$version_num""
Usage: ./STARE.sh
[-a/--annotation gene annotation file in gtf-format]
[-g/--genome input fasta file in RefSeq format]
[-s/--pscm file with PSCMs in transfac format] OR [-p/--psem file with PSEMs of TFs]
[-o/--output prefix_path of output files]\n
Optional parameters:
[-b/--bed_file bed file containing open chromatin regions, leave empty to run with promoter windows]
[-w/--window window size around TSS for mapping regions to genes (default 50KB; 5MB for ABC-mode)]
[-n/--column column in the -b file containing the average per base signal within a peak, start counting at 1]
[-c/--cores number of cores to use (default 1)]
[-x/--exclude_bed bed-file with regions to exclude (e.g. blacklisted regions)]
[-u/--genes file with rows of gene IDs/symbols to limit the output (else all in gtf)]
[-i/--tss_mode 'all_tss' to average across all annotated TSS for ABC-scoring or '5_tss' to use only the 5' TSS (default all_tss)]
[-q/--adapted_abc whether to use the adapted ABC-score (default True)]
[-f/--contact_folder folder with normalized Hi-C contact files for each chromosome in coordinate format, expects gzipped files. Set to False to use a contact estimate based on distance.]
[-k/--bin_size bin-size of the Hi-C files]
[-t/--cutoff cut-off for the ABC-score (default 0.02), set to 0 to get all scored interactions]
[-d/--pseudocount whether to use pseudocount for the contact frequency in the ABC-model (default True)]
[-m/--enhancer_window window size around enhancers for the -q adjustment (default 5MB, minimally set to -w)]
[-r/--existing_abc ABC-scoring file, if already calculated once for this input to avoid redundant calculation]
[-y/--gc_content Mean GC-content to calculate PSEMs, by default this is automatically derived from your bed_file]
[-e/--decay indicating whether exponential distance decay should be used (default TRUE, but not used in ABC-mode)]
[-z/--reshape write a binary output (default False), optional input for GAZE]"

# ------------------------------------------------------------------------------------------------------
# FETCHING INPUT
# ------------------------------------------------------------------------------------------------------
print_help=0;
print_version=0;
regions=""
genome=""
annotation=""
prefixP=""
cores=1
psems=""
pscms=""
pscm_cg=""
genes="0"
column="0"
tss_mode="all_tss"
window=""
decay="TRUE"
hic_contactfolder=""
hic_binsize=""
abc_cutoff=0.02
existing_abc="0"
exclude_regions=""
reshaping="FALSE"
pseudocount="TRUE"
adjustedABC="TRUE"
enhancer_window=5000000


die() { echo "$*" >&2; exit 2; }  # complain to STDERR and exit with error
needs_arg() { if [ -z "$OPTARG" ]; then die "No arg for --$OPT option"; fi; }  # Required to enable long options.

# Parsing command line.
while getopts hvg:b:o:c:p:s:y:u:n:i:a:w:e:f:k:t:r:x:z:d:q:m:-: OPT; do
  if [ "$OPT" = "-" ]; then   # long option: reformulate OPT and OPTARG
    OPT="${OPTARG%%=*}"       # extract long option name
    OPTARG="${OPTARG#$OPT}"   # extract long option argument (may be empty)
    OPTARG="${OPTARG#=}"      # if long option argument, remove assigning `=`
  fi
  case "$OPT" in
  h | help)	print_help=1;;
  v | version)  print_version=1;;
  g | genome)	needs_arg; genome=$OPTARG;;
	b | bed_file)	needs_arg; regions=$OPTARG;;
	o | output)	needs_arg; prefixP=$OPTARG;;
	c | cores)	needs_arg; cores=$OPTARG;;
	p | psem)	needs_arg; psems=$OPTARG;;
  s | pscm)	needs_arg; pscms=$OPTARG;;
  y | gc_content)	needs_arg; pscm_cg=$OPTARG;;
  u | genes) needs_arg; genes=$OPTARG;;
	n | column)	needs_arg; column=$OPTARG;;
  i | tss_mode) needs_arg; tss_mode=$OPTARG;;
	a | annotation)	needs_arg; annotation=$OPTARG;;
	w | window)	needs_arg; window=$OPTARG;;
	e | decay)	needs_arg; decay=$OPTARG;;
  f | contact_folder)  needs_arg; hic_contactfolder=$OPTARG;;
  k | bin_size)  needs_arg; hic_binsize=$OPTARG;;
  t | cutoff)  needs_arg; abc_cutoff=$OPTARG;;
  r | existing_abc)  needs_arg; existing_abc=$OPTARG;;
  x | exclude_bed)  needs_arg; exclude_regions=$OPTARG;;
  z | reshape)  needs_arg; reshaping=$OPTARG;;
  d | pseudocount)  needs_arg; pseudocount=$OPTARG;;
  q | adapted_abc)  needs_arg; adjustedABC=$OPTARG;;
  m | enhancer_window)  needs_arg; enhancer_window=$OPTARG;;
  ??* ) die "Illegal option --$OPT" ;;  # bad long option
  ? ) exit 2 ;;  # bad short option (error reported via getopts)
esac
done

if [ "$print_version" -eq 1 ];
then
    echo "STARE version ""$version_num"
    exit 1;
fi

if [ $OPTIND -eq 1 ] || [ "$print_help" -eq 1 ];
then
    echo -e "$help"
    exit 1;
fi

if [ -z "$genome" ] ;
then
	echo Reference genome must be specified using the -g/--genome parameter
	exit 1;
fi

if [ -z "$annotation" ] ;
then
	echo Gene annotation must be specified using the -a/--annotation parameter
	exit 1;
fi

if [ "${regions}" == "FALSE" ] || [ "${regions}" == "False" ] || [ "${regions}" == "false" ] || [ "${regions}" == "F" ] || [ "${regions}" == "0" ] ;
then
  regions="";
fi

if [ -z "$regions" ] ;
then
  column=""  # Will serve as marker for the TF_Gene_Scorer that it's the promoter-only mode.
	echo No bed-file provided, will run on promoter windows
fi

if [ -z "$prefixP" ] ;
then
	echo Prefix of output files must be specified using the -o/--output parameter
	exit 1;
fi

if [ -z "$psems" ] && [ -z "$pscms" ];
then
	echo PSEMs must be specified using the -p/--psem parameter, or PSCMs with the -s/--pscm parameter
	exit 1;
fi

if [ -n "$hic_contactfolder" ] || [ -n "$hic_binsize" ] && [[ "$existing_abc" == "0" ]];  # With an existing ABC-file, the other flags are ignored.
then  # If the hic-contactfolder is set to false we can run the ABC-scoring with contact estimate.
  if [ "${hic_contactfolder}" == "FALSE" ] || [ "${hic_contactfolder}" == "False" ] || [ "${hic_contactfolder}" == "false" ] || [ "${hic_contactfolder}" == "F" ] || [ "${hic_contactfolder}" == "0" ];
  then
    if [ "${column}" == "0" ];
    then
      echo "For the ABC-score calculation with contact estimate the column with the peak signal (-n/--column) is required."
      exit 1;
    fi
  else
    if [ -z "$hic_contactfolder" ] || [ -z "$hic_binsize" ] || [ "${column}" == "0" ];
    then
      echo "For the ABC-score calculation the column with the peak signal (-n/--column), the path to the normalized contact files (-f/--contact_folder) as well as the the bin size (-k/--bin_size) are required."
      exit 1;
    fi
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

if [ "${decay}" == "FALSE" ] || [ "${decay}" == "False" ] || [ "${decay}" == "false" ] || [ "${decay}" == "F" ] || [ "${decay}" == "0" ] ;
then
  decay="FALSE";
fi

if [ "${reshaping}" == "TRUE" ] || [ "${reshaping}" == "True" ] || [ "${reshaping}" == "true" ] || [ "${reshaping}" == "T" ] || [ "${reshaping}" == "1" ] ;
then
  reshaping="TRUE";
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
  echo "$prefixP" " Folder already exists, remove it or change the function call prefix in the -o/--output flag"
  exit 1;
fi

mkdir "${prefixP}"
base_prefix=$(basename "${prefixP}")
prefix_path=$prefixP"/"$base_prefix
working_dir=$(cd "$(dirname "$0")" && pwd -P)

metadatafile=${prefix_path}_metadata.amd.tsv
# Create metadata file.
touch "$metadatafile"
echo "[Description]" >> "$metadatafile"
echo "process	STARE ""$version_num" >> "$metadatafile"
echo -e "run_by_user\t""$USER" >> "$metadatafile"
echo -e "date\t""$d" >> "$metadatafile"
echo -e "time\t""$t" >> "$metadatafile"
echo -e "analysis_id\t""$prefix_path" >> "$metadatafile"
echo "" >> "$metadatafile"
echo "[Command]" >> "$metadatafile"
echo "STARE.sh ""$*" >> "$metadatafile"
echo "" >> "$metadatafile"
echo "[Inputs]" >> "$metadatafile"
if [ -n "$regions" ] ;
  then
    echo -e "-b region_file\t""$regions" >> "$metadatafile"
else
    echo -e "-b region_file\t""Using promoter windows" >> "$metadatafile"
fi
if [ -n "$column" ] ;
then
	echo -e "-n signal_column\t""$column" >> "$metadatafile"
fi
if [[ "$genes" != "0" ]];
then
  echo -e "-u gene set\t""$genes" >> "$metadatafile"
fi
if [ -n "$exclude_regions" ] ;
then
	echo -e "-x excluded regions\t""$exclude_regions" >> "$metadatafile"
fi
echo "" >> "$metadatafile"
echo "[References]" >> "$metadatafile"
echo -e "-g genome_reference\t""$genome" >> "$metadatafile"
if [ -n "$pscms" ];
then
  echo -e "-s pscms\t""$pscms" >> "$metadatafile"
  if [ -n "$pscm_cg" ];
  then
    echo -e "-y GC-content\t""$pscm_cg" >> "$metadatafile"
  else
    echo -e "GC-content\tautomatic from bed_file" >> "$metadatafile"
  fi
else
  echo -e "-p psems\t""$psems" >> "$metadatafile"
fi
echo -e "-a genome_annotation\t""$annotation">> "$metadatafile"

echo "" >> "$metadatafile"
echo "[Output path]" >> "$metadatafile"
echo "$prefixP" >> "$metadatafile"

echo "" >> "$metadatafile"
echo "[Parameters]" >> "$metadatafile"
echo -e "-c cores\t""$cores" >> "$metadatafile"
echo -e "-w window\t"$window >> "$metadatafile"
echo -e "-i tss_mode\t""$tss_mode" >> "$metadatafile"
if [ -z "$hic_contactfolder" ] && [[ "$existing_abc" == "0" ]];
then
  echo -e "-e decay\t""$decay" >> "$metadatafile"
fi
if [ -n "$hic_contactfolder" ] && [[ "$existing_abc" == "0" ]];
then
  echo -e "-f path with hi-c contact files\t""$hic_contactfolder" >> "$metadatafile"
  echo -e "-k bin size of hi-c contacts\t""$hic_binsize" >> "$metadatafile"
  echo -e "-t ABC-score cut-off\t""$abc_cutoff" >> "$metadatafile"
  echo -e "-d Use pseudocount for contact frequency\t""$pseudocount" >> "$metadatafile"
  echo -e "-q Use adaptedABC version\t""$adjustedABC" >> "$metadatafile"
  echo -e "-m Window size for the adaptedABC\t""$enhancer_window" >> "$metadatafile"
fi
if [[ "$existing_abc" != "0" ]];
then
  echo -e "-r existing ABC-score file that was used\t""$existing_abc" >> "$metadatafile"
fi
if [[ "$reshaping" == "TRUE" ]];
then
  echo -e "-z Reshaping to binary output\t""$reshaping" >> "$metadatafile"
fi

echo "" >> "$metadatafile"
echo "[Metrics]" >> "$metadatafile"
if [ -n "$regions" ] ;
  then
  numReg=`grep -c . "$regions"`
  echo -e "Number of provided regions\t""$numReg" >> "$metadatafile"
fi
if [ -n "$pscms" ];
then
  numMat=`grep "//" "$pscms" | wc -l`
else
  numMat=`grep ">" "$psems" | wc -l`
fi
echo -e "Number of considered psems\t""$numMat" >> "$metadatafile"

# ------------------------------------------------------------------------------------------------------
# REGION PROCESSING
# ------------------------------------------------------------------------------------------------------
echo "Preprocessing region file"

# Check for chr-prefix in the sequence file. Only true if the first line of the fasta file starts with '>chr'.
if [ "$(sed -n '/^>chr/p;q' "$genome")" ]; then
  chrPrefix="TRUE"
else
  chrPrefix="FALSE"
fi

if [ -z "$regions" ] ;
then
  "${working_dir}"/Promoter_Windows -a "${annotation}" -w "${window}" -u "${genes}" -o ${prefixP}"/"${base_prefix}
  regions=${prefixP}"/"${base_prefix}"_STARE_PromoterWindows.bed"
fi

filteredRegions=$prefix_path"_candidate_binding_regions"
sed 's/chr//g' "$regions" >  "${filteredRegions}"_Filtered_Regions.bed
sort -s -k1,1 -k2,2 -k3,3 "${filteredRegions}"_Filtered_Regions.bed | uniq > "${filteredRegions}"_sorted.bed
rm "${filteredRegions}"_Filtered_Regions.bed

# Remove regions that overlap with regions in the $exclude_regions file with the intersect -v flag.
if [ -n "$exclude_regions" ] ;
then
  sed 's/chr//g' "$exclude_regions" > "${exclude_regions}"_noPrefix.bed
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

# Generating name of the fasta file containing the overlapping regions.
openRegionSequences=${prefix_path}.OpenChromatin.fasta
# Run bedtools to get a fasta file containing the sequence data for predicted open chromatin regions contained in the bedfile.
bedtools getfasta -fi "$genome" -bed "${getFastaRegion}" -fo "$openRegionSequences"
if [ "${chrPrefix}" == "TRUE" ];
then
	rm "${getFastaRegion}"
fi

# Replace invalid characters in the fasta file with 'N', but skip the id rows (">...").
"${working_dir}"/ReplaceInvalidChars -i "$openRegionSequences" -o "${prefix_path}"_FilteredSequences.fa -d "${prefix_path}"_SeqMeta.txt
rm "$openRegionSequences"

# Transform PSCM to PSEM, if necessary.
if [ -n "$pscms" ] ;
then
  echo "Transforming transfac-PSCMs to PSEMs"
  if [ -z "$pscm_cg" ] ;
  then
    pscm_cg=`sed -n '2p' "${prefix_path}"_SeqMeta.txt`
  fi
  psems=${prefix_path}_$(basename "$pscms").PSEM
  "${working_dir}"/PSCM_to_PSEM "${pscms}" "${pscm_cg}" > "${psems}"
fi


# ------------------------------------------------------------------------------------------------------
# TF-AFFINITY WITH TRAP
# ------------------------------------------------------------------------------------------------------
startt=`date +%s`
# Use TRAP to compute transcription factor affinities to the above extracted sequences.
#affinity="/projects/triangulate/work/STARE/Hocker_scHeart/Hocker_Affinities.txt"
affinity=${prefix_path}_Affinity.txt
echo "Starting TRAP"
"${working_dir}"/TRAPmulti "$psems" "${prefix_path}"_FilteredSequences.fa "${prefix_path}"_SeqMeta.txt "$cores" > "${affinity}"
rm "${prefix_path}"_FilteredSequences.fa
rm "${prefix_path}"_SeqMeta.txt

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
  "${working_dir}"/STARE_ABCpp -b "${filteredRegions}"_sorted.bed -n "${column}" -a "${annotation}" -w "${window}" -f "${hic_contactfolder}" -k "${hic_binsize}" -t "${abc_cutoff}" -o "${abc_prefix_path}" -d "${pseudocount}" -q "${adjustedABC}" -m "${enhancer_window}" -c "${cores}" -u "${genes}" -i "${tss_mode}"
  existing_abc=${abc_prefix_path}"_ABCpp_scoredInteractions.txt.gz"
fi

# ------------------------------------------------------------------------------------------------------
# GET TF-GENE AFFINITIES
# ------------------------------------------------------------------------------------------------------
# The gene view is generated.
startg=`date +%s`
echo "Generating TF-Gene scores"
mkdir "${prefixP}""/Gene_TF_matrices"
"${working_dir}"/TF_Gene_Scorer_Reshape -a "${annotation}" -b "${filteredRegions}"_sorted.bed -n "${column}" -i "${affinity}" -o "${prefixP}"/Gene_TF_matrices/"${base_prefix}" -p "${psems}" -w ${window} -e "${decay}" -c "${cores}" -abc "${existing_abc}" -z "${reshaping}" -u "${genes}"

endg=`date +%s`
echo $((endg-startg))"s TF-Gene Scores"


# Clean-Up
rm "${affinity}"
rm "${filteredRegions}"_sorted.bed

echo "Congratulations it worked!"

