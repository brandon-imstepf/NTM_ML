#! /usr/bin/env bash

if [ "$#" -ne 10 ]; then
    echo "Invalid number of arguments provided"
else
    mkdir -p ./output_data
    python ./nn-hyp-tuning.py $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
fi
