##----------------------------------------------
# Overlaying gated populations on histocytometry 
##----------------------------------------------
include("/stor/work/Ehrlich/Users/John/projects/AKOYA/scripts/flowSOM/seg_mask_coloring_functions.jl")
## Loads in cluster_mask_threaded() and required packages 


hc_dir = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images/images_2023-08-10"
hc_dirs = readdir(hc_dir, join= true)
hc_dirs = [dir for dir in hc_dirs if contains(dir, r"DAPI_Sirpa_CD63_MerTK_A") & !contains(dir, "data|plots|knowledge_tables")]

## Livecomp01 is running DAPI_Sirpa_CD63_MerTK_A

my_colors= Dict("not_annotated" =>[199.0, 33.0, 221.0],
                "aDC2" => [209.0, 74.0, 0.0],
                "Macrophage" => [0.0, 140.0, 0.0])

for histo_dir in hc_dirs
    quant_path = joinpath(histo_dir, "combined_gated_cell_types.csv")
    mask_path  = joinpath(histo_dir, "cp2_mask.tif")
    out_path   = joinpath(histo_dir, "gating_pops.tif")
                
    cluster_mask_histocytometry(mask_path,
                                quant_path, 
                                out_path,
                                my_colors,
                                (8_000, 8_000),
                                (13_000, 13_000),
                                true)
end
