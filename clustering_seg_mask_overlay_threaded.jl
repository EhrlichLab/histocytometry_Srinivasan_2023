##----------------------------------------------
# Overlaying gated populations on histocytometry 
##----------------------------------------------
include("/stor/work/Ehrlich/Users/John/projects/AKOYA/scripts/flowSOM/seg_mask_coloring_functions.jl")
## Loads in cluster_mask_threaded() and required packages 


# hc_dir = "/stor/work/Ehrlich/Users/John/projects/misc/histocytometry"
hc_dir= "/stor/scratch/Ehrlich/Users/John/CCR4/T_cell_zone"
hc_dirs = readdir(hc_dir, join= true)
hc_dirs= [dir for dir in hc_dirs if contains(dir, "I-LN")]

# cluster_ids = ["not_annotated", "cDC2", "CD11c+_Sirpa-", "aDC2"]
# my_colors = distinguishable_colors(length(cluster_ids), 
#                                    [RGB(1,1,1), RGB(0,0,0)], 
#                                    dropseed = true)
# my_colors = channelview(my_colors) .* 255

# my_colors= Dict("not_annotated" => [199.0, 33.0, 221.0], ## Pink/Magenta
#                 "cDC2"          => [209.0, 74.0, 0.0],   ## Orangish
#                 "CD11c+_Sirpa-" => [0.0, 140.0, 0.0],    ## Green  
#                 "aDC2"          => [0.0, 127.0, 177.0])  ## Blue

## I need to make the combined_gated_cellt_types.csv or something similar for the CCR4 images. 
## THe combined_gated_cell_types.csv files has cell_type, CellID, and X_centroid, Y_centroid
my_colors= Dict("not_annotated" => [199.0, 33.0, 221.0], ## Pink/Magenta
                "CD4SP"          => [209.0, 74.0, 0.0],   ## Orangish
                "CD8SP" => [0.0, 140.0, 0.0],    ## Green  
                "Bcell"          => [0.0, 127.0, 177.0])  ## Blue
    
for histo_dir in hc_dirs
    # histo_dir  = "/stor/work/Ehrlich/Users/John/projects/misc/histocytometry/noBack_pan_test"
    # histo_dir  = "/stor/work/Ehrlich/Users/John/projects/misc/histocytometry/DAPI_CD63_whole_cell0_40x_Cortex_DAPI_CD63_CD11c_Sirpa"
    quant_path = joinpath(histo_dir, "combined_gated_cell_types.csv")
    mask_path  = joinpath(histo_dir, "mask.tif")
    out_path   = joinpath(histo_dir, "gating_pops.tif")
                
    cluster_mask_histocytometry(mask_path,
                                quant_path, 
                                out_path,
                                my_colors,
                                (2000, 2000),
                                (3000, 3000),
                                true)
end
