# Histocytometry Srinivasan et al., 2023


### Main Code Outline:

**Thymic region segmentation**

1. segment_thymic_region.ipynb: Segmentation of whole-lobe, medulla, and cortex

2. segment_and_draw_CMJ.ipynb: segmentation of CMJ and drawing CMJ for manuscript images 

3. find_thymus_region_area.ipynb: 

**Single-cell segmentation**

1. making_histocytometry_thymus_cellpose_training_data.ipynb: Manually curated training regions for custom Cellpose model

2. DC_Mac_segmentation.py: Single-cell segmentation of images using custom Cellpose model

3. per_cell_marker_quantification.sh: Extract marker expression values from segmented cells

    a. filter_cell_area_size.R

    b. CLI_csv_to_fcs.R
    

**Interpretaion**

1. combining_gated_hc_pops.R

2. Join_gated_data_w_new_CMJ_region.R

2. validate_gating_w_images.ipynb


**Manuscript figures**

1. histocytometry_paper_plots.qmd





