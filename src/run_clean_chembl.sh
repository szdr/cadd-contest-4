#!/bin/bash

input_files=../data/Chembl_bioactivity*.tsv
for input_file in $input_files; do
    output_file=`echo $input_file | sed -e "s/.tsv/_cleaned.csv/" | sed -e "s/data/output/"`
    python clean_chembl.py $input_file $output_file
done
