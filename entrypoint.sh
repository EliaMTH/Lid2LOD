#!/bin/bash
set -e

# Input arguments passed to the container
INPUT1="$1" # <building_footprints.shp>  
INPUT2="$2" # <points.las>
INPUT4="$3" # <output_folder> 

# Internal working folder (hidden from the user)
WORKING_FOLDER="/working_folder"
mkdir -p "$WORKING_FOLDER"

echo "---------------------------------------"
echo "Running pt1..."
echo "---------------------------------------"
python3 /python1/Lid2LODpt1.py "$INPUT1" "$INPUT2" "$WORKING_FOLDER"
python3 /python1_1/Lid2LODpt1b.py "$INPUT2" "$WORKING_FOLDER/ground_polygon.off"

echo "---------------------------------------"
echo "Running pt2..."
echo "---------------------------------------"
/CPP1/triangulate_city/build/triangulate_city \
    "$WORKING_FOLDER/ground_polygon.off" \
    "." \
    "$WORKING_FOLDER/buildings" \
    "$INPUT4" \
    0

echo "---------------------------------------"
echo "Running pt3..."
echo "---------------------------------------"
python3 /python2/Lid2LODpt2.py "$WORKING_FOLDER" "$INPUT4"

echo "---------------------------------------"
echo "Pipeline complete!"
echo "---------------------------------------"
