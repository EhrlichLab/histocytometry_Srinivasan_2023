##----------------------------------------------
# Overlaying gated populations on histocytometry 
##----------------------------------------------
include("/stor/work/Ehrlich/Users/John/projects/AKOYA/scripts/flowSOM/seg_mask_coloring_functions.jl")
## Loads in cluster_mask_threaded() and required packages 


hc_dir = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images/images_2023-08-10"
hc_dirs = readdir(hc_dir, join= true)
hc_dirs = [dir for dir in hc_dirs if contains(dir, r"SiglecH_A$") & !contains(dir, "data|plots|knowledge_tables")]

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
# my_colors= Dict("not_annotated" => [199.0, 33.0, 221.0], ## Pink/Magenta
#                 "CD4SP"          => [209.0, 74.0, 0.0],   ## Orangish
#                 "CD8SP"          => [0.0, 140.0, 0.0],    ## Green  
#                 "Bcell"          => [0.0, 127.0, 177.0])  ## Blue
    
my_colors= Dict("not_annotated" =>[199.0, 33.0, 221.0],
                "SiglecH+" => [209.0, 74.0, 0.0],
                "B220+" => [0.0, 140.0, 0.0],
                "pDC" => [0.0, 127.0, 177.0])

for histo_dir in hc_dirs
    quant_path = joinpath(histo_dir, "combined_gated_cell_types.csv")
    mask_path  = joinpath(histo_dir, "cp2_mask.tif")
    out_path   = joinpath(histo_dir, "gating_pops.tif")
                
    cluster_mask_histocytometry(mask_path,
                                quant_path, 
                                out_path,
                                my_colors,
                                (2_500, 17_500),
                                (6_000, 19_000),
                                true)
end
