#!/bin/bash

input_files=../output/SIRT*_random_sampled_10.csv
for input_file in $input_files; do
    echo $input_file
    output_file=`echo $input_file | sed -e "s/.csv/.svmlight/"`
    python make_feature_smi.py $input_file $output_file
done
