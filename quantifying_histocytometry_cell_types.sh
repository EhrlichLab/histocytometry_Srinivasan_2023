#!/bin/bash

## Loop through histocytometry images
## Find cell type masks for each image

hc_dir=/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images
cd $hc_dir

image_dirs=$(ls | grep "0x" | grep "pan")
## We're only using panoramas for the final analyses, so we can remove the rest. 

for image_dir in $image_dirs; do
    cell_type_dirs=$(ls $image_dir | grep "^cell_type_")
    for cell_type_dir in $cell_type_dirs; do
        ## Running pipeline
        /stor/work/Ehrlich/Users/John/projects/misc/histocytometry/scripts/cell_type_quantification_pipeline.sh $hc_dir/$image_dir/$cell_type_dir
    done
    ## Combine cell type files 
    Rscript \
    /stor/work/Ehrlich/Users/John/projects/misc/histocytometry/scripts/combine_cell_type_csvs.R \
    --path $image_dir
done

## Overlaying cell types on images

