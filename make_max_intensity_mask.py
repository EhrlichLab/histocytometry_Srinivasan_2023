import tifffile
import numpy as np
import os
import pandas as pd
## Using Mesmer_pypi conda environment

def make_max_intensity_mask(raw_path, out_path, mem_marker_range):
    img = tifffile.imread(raw_path)
    mem_marker_max = np.max(img[mem_marker_range,...], axis= 0)
    mem_marker_max = mem_marker_max[np.newaxis, ...] 

    max_mask = np.concatenate((img[0,np.newaxis, ...], mem_marker_max))

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    tifffile.imwrite(os.path.join(out_path, "image.tif"), max_mask)

    markers_df = pd.DataFrame(data= {
        "cycle"        : [1,1],
        "marker_name"  : ["DAPI", "max"],
        "row_num"      : [0, 1]})
    markers_df.to_csv(os.path.join(out_path, "markers.csv"), index= False)



if __name__ == "__main__":
    img_path        = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images"
    img_dirs= os.listdir(img_path)
    marker_positions= [1,2,3]
        ## This assumes that DAPI is the first marker, which it is for this histocytometry data set.

    ## Making Mask
    for img_dir in img_dirs:
        make_max_intensity_mask(raw_path = os.path.join(img_path, img_dir, "image.tif"),
                                out_path = os.path.join(img_path, img_dir, "max_membrane_mask"),
                                mem_marker_range = marker_positions)
        print("Done with one image!")

