#!/bin/bash

DIRECTORY="../Results"

mkdir $DIRECTORY
cp *.json $DIRECTORY
cp ./utils/convert_json.py $DIRECTORY
cp ./utils/plot_results.py $DIRECTORY
cd $DIRECTORY

python3 convert_json.py gr.json
python3 convert_json.py pot.json
python3 convert_json.py err.json


if [ "$1" == "skip" ]; then
    python3 plot_results.py err
    python3 plot_results.py final
else
    python3 plot_results.py all
fi

cd ..