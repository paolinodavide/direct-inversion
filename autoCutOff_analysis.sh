#!/bin/bash

# Path to the config.json file
CONFIG_FILE="config.json"

# Iterate over x_cut values from 1.0 to 2.5 spaced by 0.1
for x_cut in $(seq 2.0 0.1 4); do
    if [ -f "$CONFIG_FILE" ]; then
        # Update the 'x_cut' parameter in the config.json file
        jq --argjson x_cut "$x_cut" '.x_cut = $x_cut' "$CONFIG_FILE" > tmp.$$.json && mv tmp.$$.json "$CONFIG_FILE"
        echo "Updated 'x_cut' to $x_cut in $CONFIG_FILE"
    else
        echo "Error: $CONFIG_FILE not found!"
        exit 1
    fi

    # Run the Python script
    python3 grinter_parallel.py

    # Ask if rerun is needed or generate results and go to next
    read -p "Do you want to rerun the Python script for x_cut=$x_cut? (y/n): " rerun
    while [[ "$rerun" == "y" || "$rerun" == "Y" ]]; do
        python3 grinter_parallel.py
        read -p "Do you want to rerun the Python script for x_cut=$x_cut again? (y/n): " rerun
    done

    # Run the Bash script
    bash Cutoff_res.sh
done