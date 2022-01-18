#!/bin/bash

# Checks whether all files in a folder are identical to files with the same name in the other folder. And checks
# whether both have the same number of files.

help="[-a first folder, all files within will be compared to those with identical name in folder -b]\n
      [-b second folder]\n
      [-x files containing that string are excluded from comparison]\n"

while getopts "a:b:x:" o;
do
case $o in
	a)	folder1=$OPTARG;;
	b)	folder2=$OPTARG;;
	x)	exclude=$OPTARG;;
  *)  exit 1;;
esac
done

if [ $OPTIND -eq 1 ] ;
then
    echo -e "$help"
    exit 1;
fi

for file in $(find "$folder1" -type f); do
  if [[ $file != *.DS_Store ]] && [[ $file != *"$exclude"* ]] && [[ $file != */._* ]];
  then
    base_file=${file#"$folder1"}  # Basename doesn't work for subdirectories.
    if [[ $file == *.gz ]];
    then  # Tried zcmp, zdiff and cmp -i but neither did really work on both Mac and Ubuntu.
      gunzip "$folder1""$base_file"
      gunzip "$folder2""$base_file"
      unzipped_base=${base_file%.gz}
      sort "$folder1""$unzipped_base" > "$folder1""$unzipped_base"_sorted  # We need to sort, as we write from sets.
      sort "$folder2""$unzipped_base" > "$folder2""$unzipped_base"_sorted
      if ! cmp --silent -- "$folder1""$unzipped_base"_sorted "$folder2""$unzipped_base"_sorted; then
        echo "ERROR: test file differs with control file:" "$folder1""$base_file" "$folder2""$base_file"
      fi
      gzip "$folder1""$unzipped_base"
      gzip "$folder2""$unzipped_base"
      rm  "$folder1""$unzipped_base"_sorted
      rm  "$folder2""$unzipped_base"_sorted
    else
      sort "$folder1""$base_file" > "$folder1""$base_file"_sorted
      sort "$folder2""$base_file" > "$folder2""$base_file"_sorted
      if ! cmp --silent -- "$folder1""$base_file" "$folder2""$base_file"; then
        echo "ERROR: test file differs with control file:" "$folder1""$base_file" "$folder2""$base_file"
      fi
      rm  "$folder1""$base_file"_sorted
      rm  "$folder2""$base_file"_sorted
    fi
  fi
done


folder1_num=$(ls -1q "$folder1"/* | grep -v "^$" | wc -l)
folder2_num=$(ls -1q "$folder2"/* | grep -v "^$" | wc -l)

if [[ $folder1_num != $folder2_num ]];
then
  echo "ERROR: the number of files does not agree:" "$folder1" "$folder2"
fi
