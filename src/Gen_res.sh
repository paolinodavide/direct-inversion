#!/bin/bash

DIRECTORY="Results"

mkdir $DIRECTORY
cp *.json $DIRECTORY
cp convert_json.py $DIRECTORY
cd $DIRECTORY

python3 convert_json.py gr.json
python3 convert_json.py pot.json
python3 convert_json.py err.json

cd ..
