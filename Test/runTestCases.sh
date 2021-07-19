echo TestV1: Windows 3kb - Decay
bash ../Code/STARE.sh -g example_sequence.fa -b example_regions.bed  -o Test_V1 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000
echo ""

echo TestV2: Windows 3kb - No Decay
bash ../Code/STARE.sh -g example_sequence.fa -b example_regions.bed  -o Test_V2 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE
echo ""

echo TestV3: Windows 3kb - No Decay - Signal Feature
bash ../Code/STARE.sh -g example_sequence.fa -b example_regions.bed  -o Test_V3 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -n 4
echo ""

echo TestV4: Windows 3kb - Decay - Chr prefix in the reference genome
bash ../Code/STARE.sh -g example_sequence_chr.fa -b example_regions.bed  -o Test_V4 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -j
echo ""

echo TestV5: Windows 3kb - Decay - Exclude regions
bash ../Code/STARE.sh -g example_sequence.fa -b example_regions.bed  -o Test_V5 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -x example_exclude_regions.bed
echo ""

echo TestV6: Windows 3kb - No Decay - Multiple Signal Columns
bash ../Code/STARE.sh -g example_sequence.fa -b example_regions.bed  -o Test_V6 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a example_annotation.gtf  -w 3000 -e FALSE -n 4+
echo ""

echo TestV7: ABC Windows 500kb - Multiple Signal Columns - ABC-annotation - ChrPrefix - Reshaping
bash ../Code/STARE.sh -g example_sequence_chr21part.fa -b ABC_example_regions.bed -o Test_V7 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a ABC_example_annotation.gtf -w 500000 -e False -n 4+ -t 0 -j -f ABC_example_HiCFiles/ -k 5000 -t 0 -z
echo ""

echo TestV8: Existing ABC-annotation - ChrPrefix
bash ../Code/STARE.sh -g example_sequence_chr21part.fa -b ABC_example_regions.bed  -o Test_V8 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a ABC_example_annotation.gtf -e False -j -r example_ABCpp_scoredInteractions_c4.txt.gz
echo ""

echo TestV9: Multiple Existing ABC-annotations - ChrPrefix
bash ../Code/STARE.sh -g example_sequence_chr21part.fa -b ABC_example_regions.bed  -o Test_V9 -p ../PWMs/2.0/human_jaspar_hoc_kellis.PSEM  -a ABC_example_annotation.gtf -e False -j -r example_ABCpp_scoredInteractions_c4.txt.gz -n 4-5
echo ""

echo TestV10: separate ABC scoring
../Code/STARE_ABCpp_absolute -b ABC_example_regions_ChrPrefix.bed -n 4 -a ABC_example_annotation.gtf -gw 500000 -bin 5000 -t 0 -cf ABC_example_HiCFiles/ -d TestV10
echo ""

echo TestV11: separate ABC scoring - Multiple Signal Columns
../Code/STARE_ABCpp_absolute -b ABC_example_regions.bed -n 4-5 -a ABC_example_annotation.gtf -gw 500000 -bin 5000 -t 0 -cf ABC_example_HiCFiles/ -d TestV11
echo ""
