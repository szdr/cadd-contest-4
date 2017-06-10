#!/bin/bash

input_files=../output/Chembl_bioactivity_SIRT*_cleaned.svmlight
sample_ratio=10
for input_file in $input_files; do
    echo $input_file
    output_file=`echo $input_file | awk -F'_' '{print $3}'`_random_sampled_$sample_ratio.csv
    python random_sample_from_zinc.py $input_file $output_file --sample_ratio $sample_ratio
done
