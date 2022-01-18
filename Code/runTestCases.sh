working_dir=$(cd "$(dirname "$0")" && pwd -P)
test_dir="$working_dir"/../Test/
pwm_file="$working_dir"/../PWMs/2.0/human_jaspar_hoc_kellis.PSEM


echo TestV1: Windows 3kb - Decay
bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence.fa -b "$test_dir"/Test_Data/example_regions.bed  -o "$test_dir"/Test_V1 -p "$pwm_file" -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V1 -b "$test_dir"/Test_Controls/Test_V1 -x metadata
echo ""

echo TestV2: Windows 3kb - No Decay
bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence.fa -b "$test_dir"/Test_Data/example_regions.bed  -o "$test_dir"/Test_V2 -p "$pwm_file" -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000 -e FALSE
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V2 -b "$test_dir"/Test_Controls/Test_V2 -x metadata
echo ""

echo TestV3: Windows 3kb - No Decay - Signal Feature
bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence.fa -b "$test_dir"/Test_Data/example_regions.bed  -o "$test_dir"/Test_V3 -p "$pwm_file" -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000 -e FALSE -n 4
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V3 -b "$test_dir"/Test_Controls/Test_V3 -x metadata
echo ""

echo TestV4: Windows 3kb - Decay - Chr prefix in the reference genome
bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr.fa -b "$test_dir"/Test_Data/example_regions.bed  -o "$test_dir"/Test_V4 -p "$pwm_file" -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000 -j True
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V4 -b "$test_dir"/Test_Controls/Test_V4 -x metadata
echo ""

echo TestV5: Windows 3kb - Decay - Exclude regions
bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence.fa -b "$test_dir"/Test_Data/example_regions.bed  -o "$test_dir"/Test_V5 -p "$pwm_file" -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000 -x "$test_dir"/Test_Data/example_exclude_regions.bed
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V5 -b "$test_dir"/Test_Controls/Test_V5 -x metadata
echo ""

echo TestV6: Windows 3kb - No Decay - Multiple Signal Columns
bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence.fa -b "$test_dir"/Test_Data/example_regions.bed  -o "$test_dir"/Test_V6 -p "$pwm_file" -a "$test_dir"/Test_Data/example_annotation.gtf -w 3000 -e FALSE -n 4+
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V6 -b "$test_dir"/Test_Controls/Test_V6 -x metadata
echo ""

echo TestV7: ABC Windows 500kb - Multiple Signal Columns - ABC-annotation - ChrPrefix
bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr21part.fa -b "$test_dir"/Test_Data/ABC_example_regions.bed -o "$test_dir"/Test_V7 -p "$pwm_file" -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -w 500000 -n 4+ -t 0 -j True -f "$test_dir"/Test_Data/ABC_example_HiCFiles/ -k 5000 -t 0
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V7 -b "$test_dir"/Test_Controls/Test_V7 -x metadata
echo ""

echo TestV8: Existing ABC-annotation - ChrPrefix
bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr21part.fa -b "$test_dir"/Test_Data/ABC_example_regions.bed -n 4+ -o "$test_dir"/Test_V8 -p "$pwm_file" -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -j True -r "$test_dir"/Test_Data/example_ABCpp_scoredInteractions_c4.txt.gz
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V8 -b "$test_dir"/Test_Controls/Test_V8 -x metadata
echo ""

echo TestV9: ABC Windows 500kb NOT adapted - Multiple Signal Columns - ABC-annotation - ChrPrefix
bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr21part.fa -b "$test_dir"/Test_Data/ABC_example_regions.bed -o "$test_dir"/Test_V9 -p "$pwm_file" -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -w 500000 -n 4+ -t 0 -j True -f "$test_dir"/Test_Data/ABC_example_HiCFiles/ -k 5000 -t 0 -q False
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V9 -b "$test_dir"/Test_Controls/Test_V9 -x metadata
echo ""

echo TestV10: Multiple Existing ABC-annotations - ChrPrefix
bash "$working_dir"/STARE.sh -g "$test_dir"/Test_Data/example_sequence_chr21part.fa -b "$test_dir"/Test_Data/ABC_example_regions.bed  -o "$test_dir"/Test_V10 -p "$pwm_file" -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -e False -j True -r "$test_dir"/Test_Data/example_ABCpp_scoredInteractions_c4.txt.gz -n 4-5
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V10 -b "$test_dir"/Test_Controls/Test_V10 -x metadata
echo ""

echo TestV11: separate ABC scoring
mkdir "$test_dir"/Test_V11
"$working_dir"/STARE_ABCpp -b "$test_dir"/Test_Data/ABC_example_regions_ChrPrefix.bed -n 4 -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -d "$test_dir"/Test_V11/Test_V11 -gw 500000 -bin 5000 -t 0 -cf "$test_dir"/Test_Data/ABC_example_HiCFiles/
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V11 -b "$test_dir"/Test_Controls/Test_V11 -x metadata
echo ""

echo TestV12: separate ABC scoring - Multiple Signal Columns
mkdir "$test_dir"/Test_V12
"$working_dir"/STARE_ABCpp -b "$test_dir"/Test_Data/ABC_example_regions.bed -n 4-5 -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -d "$test_dir"/Test_V12/Test_V12 -gw 500000 -bin 5000 -t 0 -cf "$test_dir"/Test_Data/ABC_example_HiCFiles/
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V12 -b "$test_dir"/Test_Controls/Test_V12 -x metadata
echo ""

echo TestV13: separate ABC scoring NOT adapted - Multiple Signal Columns
mkdir "$test_dir"/Test_V13
"$working_dir"/STARE_ABCpp -b "$test_dir"/Test_Data/ABC_example_regions.bed -n 4-5 -a "$test_dir"/Test_Data/ABC_example_annotation.gtf -d "$test_dir"/Test_V13/Test_V13 -gw 500000 -bin 5000 -t 0 -cf "$test_dir"/Test_Data/ABC_example_HiCFiles/ -q False
bash "$working_dir"/controlTestCases.sh -a "$test_dir"/Test_V13 -b "$test_dir"/Test_Controls/Test_V13 -x metadata
echo ""




