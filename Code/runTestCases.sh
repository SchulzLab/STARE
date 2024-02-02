#!/bin/bash

help="
[-i Test folder containing the subdirectories Test_Data and Test_Controls]\n
[-o directory to write the output to, each test will get its own folder there]\n
[-t leave empty or write 'all' to run all test cases, or give a number of individual ones]"

print_help=0
while getopts "hi:o:t:" o;
do
case $o in
  h)	print_help=1;;
	i)	test_dir=$OPTARG;;
	o)	output=$OPTARG;;
  t)  which_test=$OPTARG;;
  *)  exit 1;;
esac
done

if [ $OPTIND -eq 1 ] || [ "$print_help" -eq 1 ];
then
    echo -e "$help"
    exit 1;
fi


if [ -z "$test_dir" ];
then
  echo "Missing -i test folder"
  exit 1;
fi

if [ -z "$output" ];
then
  echo "Missing -o output path"
  exit 1;
fi

if [ -z "$which_test" ];
then
  which_test="all"
fi

working_dir=$(cd "$(dirname "$0")" && pwd -P)

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "1" ]];
then
  echo TestV1: Windows 3kb - Decay
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence.fa -b "$test_dir"/Test_Data/example_regions.bed -o "$output"/Test_V1 -p "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_PSEMs.txt -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V1 -b "$test_dir"/Test_Controls/Test_V1 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "2" ]];
then
  echo TestV2: Windows 3kb - No Decay
  bash "$working_dir"/STARE.sh --genome="$test_dir"/Test_Data/example_sequence.fa --bed_file="$test_dir"/Test_Data/example_regions.bed --output="$output"/Test_V2 --psem="$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_PSEMs.txt --annotation="$test_dir"/Test_Data/example_annotation.gtf --window=3000 --decay=FALSE
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V2 -b "$test_dir"/Test_Controls/Test_V2 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "3" ]];
then
  echo TestV3: Windows 3kb - No Decay - Signal Feature
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence.fa -b "$test_dir"/Test_Data/example_regions.bed -o "$output"/Test_V3 -p "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_PSEMs.txt -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000 -e FALSE -n 4
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V3 -b "$test_dir"/Test_Controls/Test_V3 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "4" ]];
then
  echo TestV4: Windows 3kb - Sequence file has chr-prefix
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr.fa -b "$test_dir"/Test_Data/example_regions.bed -o "$output"/Test_V4 -p "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_PSEMs.txt -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V4 -b "$test_dir"/Test_Controls/Test_V4 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "5" ]];
then
  echo TestV5: Windows 3kb - Decay - Exclude regions
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence.fa -b "$test_dir"/Test_Data/example_regions.bed -o "$output"/Test_V5 -p "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_PSEMs.txt -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000 -x "$test_dir"/Test_Data/example_exclude_regions.bed
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V5 -b "$test_dir"/Test_Controls/Test_V5 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "6" ]];
then
  echo TestV6: Windows 3kb - No Decay - Multiple Signal Columns
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence.fa -b "$test_dir"/Test_Data/example_regions.bed -o "$output"/Test_V6 -p "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_PSEMs.txt -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000 -e FALSE -n 4+
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V6 -b "$test_dir"/Test_Controls/Test_V6 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "7" ]];
then
  echo TestV7: ABC 500kb - Multiple Signal Columns
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr21part.fa -b "$test_dir"/Test_Data/ABC_example_regions.bed -o "$output"/Test_V7 -p "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_PSEMs.txt -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -w 500000 -n 4+ -f "$test_dir"/Test_Data/ABC_example_HiCFiles/ -k 5000 -t 0 -i 5_tss
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V7 -b "$test_dir"/Test_Controls/Test_V7 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "8" ]];
then
  echo TestV8: ABC 500kb - Multiple Signal Columns - Gene subset
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr21part.fa -b "$test_dir"/Test_Data/ABC_example_regions.bed -o "$output"/Test_V8 -p "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_PSEMs.txt -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -w 500000 -n 4+ -f "$test_dir"/Test_Data/ABC_example_HiCFiles/ -k 5000 -t 0 -u "$test_dir"/Test_Data/example_genes.txt -i 5_tss
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V8 -b "$test_dir"/Test_Controls/Test_V8 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "9" ]];
then
  echo TestV9: Existing ABC-annotation
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr21part.fa -b "$test_dir"/Test_Data/ABC_example_regions.bed -o "$output"/Test_V9 -p "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_PSEMs.txt -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -n 4+ -r "$test_dir"/Test_Data/example_ABCpp_scoredInteractions_c4.txt.gz
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V9 -b "$test_dir"/Test_Controls/Test_V9 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "10" ]];
then
  echo TestV10: ABC 500kb NOT adapted - Multiple Signal Columns
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr21part.fa -b "$test_dir"/Test_Data/ABC_example_regions.bed -o "$output"/Test_V10 -p "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_PSEMs.txt -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -w 500000 -n 4+ -t 0 -f "$test_dir"/Test_Data/ABC_example_HiCFiles/ -k 5000 -t 0 -q False -i 5_tss
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V10 -b "$test_dir"/Test_Controls/Test_V10 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "11" ]];
then
  echo TestV11: Multiple Existing ABC-annotations
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr21part.fa -b "$test_dir"/Test_Data/ABC_example_regions.bed -o "$output"/Test_V11 -p "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_PSEMs.txt -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -e False -r "$test_dir"/Test_Data/example_ABCpp_scoredInteractions_c4.txt.gz -n 4-5 -i 5_tss
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V11 -b "$test_dir"/Test_Controls/Test_V11 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "12" ]];
then
  echo TestV12: separate ABC scoring
  mkdir "$output"/Test_V12
  "$working_dir"/STARE_ABCpp -b "$test_dir"/Test_Data/ABC_example_regions_ChrPrefix.bed -n 4 -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -o "$output"/Test_V12/Test_V12 -w 500000 -k 5000 -t 0 -f "$test_dir"/Test_Data/ABC_example_HiCFiles/ -i 5_tss
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V12 -b "$test_dir"/Test_Controls/Test_V12 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "13" ]];
then
  echo TestV13: separate ABC scoring - Multiple Signal Columns
  mkdir "$output"/Test_V13
  "$working_dir"/STARE_ABCpp -b "$test_dir"/Test_Data/ABC_example_regions.bed -n 4-5 -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -o "$output"/Test_V13/Test_V13 -w 500000 -k 5000 -t 0 -f "$test_dir"/Test_Data/ABC_example_HiCFiles/ -i 5_tss
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V13 -b "$test_dir"/Test_Controls/Test_V13 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "14" ]];
then
  echo TestV14: separate ABC scoring NOT adapted - Multiple Signal Columns
  mkdir "$output"/Test_V14
  "$working_dir"/STARE_ABCpp -b "$test_dir"/Test_Data/ABC_example_regions.bed -n 4,5 -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -o "$output"/Test_V14/Test_V14 -w 500000 -k 5000 -t 0 -f "$test_dir"/Test_Data/ABC_example_HiCFiles/ -q False -i 5_tss
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V14 -b "$test_dir"/Test_Controls/Test_V14 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "15" ]];
then
  echo TestV15: Windows 3kb - Decay - transfac PSCM file
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence.fa -b "$test_dir"/Test_Data/example_regions.bed -o "$output"/Test_V15 -s "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_transfac.txt -y 0.41 -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000 -x "$test_dir"/Test_Data/example_exclude_regions.bed -i 5_tss
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V15 -b "$test_dir"/Test_Controls/Test_V15 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "16" ]];
then
  echo TestV16: ABC 500kb - Multiple Signal Columns - transfac PSCM file with automatic base content  - Binary output
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr21part.fa -b "$test_dir"/Test_Data/ABC_example_regions.bed -o "$output"/Test_V16 -s "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_transfac.txt -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -w 500000 -n 4+ -f "$test_dir"/Test_Data/ABC_example_HiCFiles/ -k 5000 -t 0 -z True -i 5_tss
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "17" ]];
then
  echo TestV17: separate ABC scoring - Multiple Signal Columns - average TSS
  mkdir "$output"/Test_V17
  "$working_dir"/STARE_ABCpp -b "$test_dir"/Test_Data/ABC_example_regions.bed -n 4-5 -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -o "$output"/Test_V17/Test_V17 -w 500000 -k 5000 -t 0 -f "$test_dir"/Test_Data/ABC_example_HiCFiles/ -i all_tss
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V17 -b "$test_dir"/Test_Controls/Test_V17 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "18" ]];
then
  echo TestV18: separate ABC scoring - Use contact estimate - average TSS
  mkdir "$output"/Test_V18
  "$working_dir"/STARE_ABCpp -b "$test_dir"/Test_Data/ABC_example_regions.bed -n 4 -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -o "$output"/Test_V18/Test_V18 -w 500000 -t 0 -f False -i all_tss
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V18 -b "$test_dir"/Test_Controls/Test_V18 -x metadata
  echo ""
fi

if [[ "$which_test" == "all" ]] || [[ "$which_test" == "19" ]];
then
  echo TestV19: Windows 3kb - Promoter windows instead of bed-file
  bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence.fa -o "$output"/Test_V19 -s "$test_dir"/Test_Data/Jaspar_Hocomoco_Kellis_human_transfac.txt -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000 -x "$test_dir"/Test_Data/example_exclude_regions.bed
  bash "$working_dir"/controlTestCases.sh -a "$output"/Test_V19 -b "$test_dir"/Test_Controls/Test_V19 -x metadata
  echo ""
fi
