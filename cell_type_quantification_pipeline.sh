#!/bin/bash

##--------------
# Script Outline
##--------------
#1. Take cell type mask and pass .ome.tif to mcquant
#2. Combine extracted cell type data
#3. Filter out "cells" that are too small to be real. 

## Requirements: 
    ## A path to a directory that has one .tif and and one markers.csv file 
    ## I'll write csv and fcs to that directory. 


dir=$1
cd $dir
echo $dir


OME_TIFF_name=$(ls .. | grep "sub_image.ome.tif")
mask_name=$(ls | grep "_mask.tif")
marker_file=$(ls .. | grep "markers.csv")

##---------------
# Running MCquant
##---------------
/stor/home/jfm2773/anaconda3/envs/mcquant/bin/python \
    /stor/work/Ehrlich/Users/John/projects/AKOYA/scripts/cell_segmentation/mcquant/CommandSingleCellExtraction.py \
    --masks $dir/$mask_name \
    --image $dir/../$OME_TIFF_name \
    --output $dir \
    --channel_names $dir/../$marker_file \
    --intensity_props mean_intensity 

echo "Done with MCquant"


##----------------------
# Filter out small cells
##----------------------
quant_csv=$(ls | grep _mask.csv)
min_area=150

Rscript \
    /stor/work/Ehrlich/Users/John/projects/misc/histocytometry/scripts/filter_cell_area_size.R \
    --path $dir/$quant_csv \
    --threshold $min_area

echo "Done filtering out cell fragments"


##------------------
# Convert CSV to FCS
##------------------
filtered_csv=$(ls | grep "_mask_area_filtered.csv")

Rscript \ 
    /stor/work/Ehrlich/Users/John/projects/AKOYA/scripts/CLI_csv_to_fcs.R \
    --path $dir/$filtered_csv

echo "Done writing FCS"


##-------------
# Example Usage
##-------------




