#!/bin/bash

DIRECTORY="../Results_cutoff"

mkdir $DIRECTORY
cp *.json $DIRECTORY
cp ./utils/convert_json.py $DIRECTORY
cp ./utils/create_cutoff_res.py $DIRECTORY
cd $DIRECTORY

python3 convert_json.py gr.json
python3 convert_json.py pot.json

python3 create_cutoff_res.py 

find . -type f ! -name 'cutoff_results.dat' -delete