#!/usr/bin/bash


raw_dir="/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images/images_2023-08-10"
cd $raw_dir

img_dirs=$(ls ${raw_dir} | grep -E "_[A-D]$") 
mask_name="cp2_mask.tif"
marker_file="markers_tissue_50um.csv"
OME_TIFF_name="reordered_image_w_thymus_regions_50um.ome.tif"


for img_dir in $img_dirs; do
    
    ## Running MCquant w/ mcquant conda environment
    /stor/home/jfm2773/anaconda3/envs/mcquant/bin/python \
        /stor/work/Ehrlich/Users/John/projects/AKOYA/scripts/cell_segmentation/mcquant/CommandSingleCellExtraction.py \
        --masks $img_dir/$mask_name \
        --image $img_dir/$OME_TIFF_name \
        --output $img_dir \
        --channel_names $img_dir/$marker_file \
        --intensity_props intensity_median 
            ## mean_intensity works as well 

    ## Filter out small non-cellular blobs
    min_area=150
    quant_csv="reordered_image_w_thymus_regions_50um_cp2_mask.csv"
    /stor/home/jfm2773/anaconda3/envs/R_flow/bin/Rscript \
        /stor/work/Ehrlich/Users/John/projects/misc/histocytometry/scripts/filter_cell_area_size.R \
        --path $raw_dir/$img_dir/$quant_csv \
        --threshold $min_area

    # Write and scale .fcs files 
    filtered_csv="reordered_image_w_thymus_regions_50um_cp2_mask_min_area_filtered.csv"
    /stor/home/jfm2773/anaconda3/envs/R_flow/bin/Rscript \ 
        /stor/work/Ehrlich/Users/John/projects/AKOYA/scripts/CLI_csv_to_fcs.R \
        --path $raw_dir/$img_dir/$filtered_csv

    echo "Done writing FCS"

done
