##--------------------------------------------
# Running thymus histocytometry image overlays
##--------------------------------------------

include("/stor/work/Ehrlich/Users/John/projects/AKOYA/scripts/flowSOM/seg_mask_coloring_functions.jl")

raw_dir  = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images"
hc_dirs = readdir(raw_dir, join= true)
hc_dirs = [dir for dir in hc_dirs if occursin(r"(2|4)0x_", dir) & !occursin(r"pan|SOCS2|40x_pan_DAPI_CD63_CD11c_XCR1", dir) & occursin(r"DAPI_CD63_CD11c_XCR1|DAPI_CD63_CD11c_Sirpa", dir)]
    ## Checking that it is an image directory that we want to analyze (not SOCS2)
    ## 40x_pan_DAPI_CD63_CD11c_XCR1 doesn't have a segmentation mask for some reason. I'm not sure why. 

for hc_dir in hc_dirs
    mask_path    = joinpath(hc_dir, "mask.tif")
    quant_path   = joinpath(hc_dir, "scyan_pops.csv")
    out_path     = joinpath(hc_dir, "scyan_pops.tif")

    start_point= (1,1)
    end_point= size(rawview(channelview(load(mask_path))))
    cluster_mask_cell_type_threaded(quant_path, mask_path, out_path, start_point, end_point)
    println("Done with " * hc_dir * "!")
end



