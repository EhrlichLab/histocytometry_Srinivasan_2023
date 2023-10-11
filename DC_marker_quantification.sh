#!/usr/bin/bash


raw_dir="/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images/images_2023-08-10"
cd $raw_dir

img_dirs=$(ls ${raw_dir} | grep -E "_[A-D]$" | grep -Ev "CD207_CD11c_XCR1_C") 
mask_name="cp2_mask.tif"
marker_file="markers_tissue_50um.csv"
OME_TIFF_name="reordered_image_w_thymus_regions_50um.ome.tif"

img_dirs="20x_pan_DAPI_B220_CD11c_SiglecH_A"


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
        ## Run w/ the R_flow conda environment
    min_area=150
        ## I picked this number out of a hat, but sqrt(150 pixels * 0.479 um/pixel) == 8.47 um per side of a cell (assuming a square), which is smaller than expected for a DC.
        ## There's also a ton of small cell fragment noise.
    # quant_csv=$(ls ${img_dir} | grep "_mask.csv")
    quant_csv="reordered_image_w_thymus_regions_50um_cp2_mask.csv"
        ## This is restrictive.  
    /stor/home/jfm2773/anaconda3/envs/R_flow/bin/Rscript \
        /stor/work/Ehrlich/Users/John/projects/misc/histocytometry/scripts/filter_cell_area_size.R \
        --path $raw_dir/$img_dir/$quant_csv \
        --threshold $min_area

     # Write and scale .fcs files 
    # filtered_csv=$(ls ${img_dir} | grep "_min_area_filtered.csv")
    filtered_csv="reordered_image_w_thymus_regions_50um_cp2_mask_min_area_filtered.csv"
    /stor/home/jfm2773/anaconda3/envs/R_flow/bin/Rscript \ 
        /stor/work/Ehrlich/Users/John/projects/AKOYA/scripts/CLI_csv_to_fcs.R \
        --path $raw_dir/$img_dir/$filtered_csv

    echo "Done writing FCS"

    ## I need to filter by area.
done
