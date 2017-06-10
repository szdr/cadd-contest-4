#!/bin/bash

input_files=../output/Chembl_bioactivity*_cleaned.csv
for input_file in $input_files; do
    echo $input_file
    output_file=`echo $input_file | sed -e "s/.csv/.svmlight/"`
    python make_feature_chembl.py $input_file $output_file &
done
wait
