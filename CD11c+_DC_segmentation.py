##--------------------------------------------------------------------
# Running histocytometry of DC images from Cellpose 2.0 using a script
##--------------------------------------------------------------------
    ## Use cellpose2 conda environment

import numpy as np
import time, os, sys
from urllib.parse import urlparse
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
import cellpose 
from cellpose import utils, io, models
import re
import tifffile
import cv2
# from deco import concurrent, synchronized

##---------
# Functions
##---------
# @concurrent
def segment_DCs_w_cellpose(img_dir, 
                           user_model, 
                           cellprob_threshold,
                           flow_threshold,
                           img_name               = "reordered_image.ome.tif", 
                           cellpose2_channels     = [2,1], 
                           DAPI_membrane_channels = [0,2],
                           out_name               = "cp2_mask",
                           downsample             = 0
                           ) -> None:
    img     = io.imread(os.path.join(img_dir, img_name))
    sub_img = img[DAPI_membrane_channels,...]    
        ## Subsetting image to get the membrane and DAPI channel
        ## DAPI corresponds to red (1),
        ## CD11c corresponds to green (2) 

    if downsample != 0:
        sub_img = sub_img[:, 0:downsample, 0:downsample]
        img = img[:, 0:downsample, 0:downsample]

    cellpose2_channels = [2,1]
    masks, flows, styles = user_model.eval(
        x                  = sub_img, 
        channels           = cellpose2_channels, 
        cellprob_threshold = cellprob_threshold,
        flow_threshold     = flow_threshold
    )
    print("Done with modeling!")

    tifffile.imwrite(os.path.join(img_dir, out_name + ".tif"), masks)

    ## Convert mask to binary for overlaying 
    masks[masks > 0] = 1
    masks.astype("uint8")
    masks = masks[np.newaxis, ...]
    img_w_cp2_masks = np.concatenate((img, masks))
    tifffile.imwrite(os.path.join(img_dir, "reordered_image_w_" + out_name + ".ome.tif"), img_w_cp2_masks)
    print("Done saving images!")

    return None


if __name__ == "__main__":
    ## Load CD11c+ images
    raw_dir = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images/images_2023-08-30"
    CD11c_medulla_dirs = [os.path.join(raw_dir, img_dir) for img_dir in os.listdir(raw_dir) if re.search("_[A-D]$", img_dir) and not re.search("Sirpa_C$", img_dir) and re.search("medulla", img_dir)]
        ## and not re.search("MerTK", img_dir)
        ## I'm throwing in MerTK to see how they look.
    CD11c_medulla_dirs.sort()

    ## Load-in user trained model
    model_path = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images/images_2023-08-10/CP_v3_DAPI_CD11c_training/models/CP_DAPI_CD11c_v6"
    user_model = models.CellposeModel(pretrained_model= model_path)

    ## Segment DCs 
    for img_dir in CD11c_medulla_dirs:
        print(img_dir)
        segment_DCs_w_cellpose(
            img_dir            = img_dir, 
            user_model         = user_model,
            img_name           = "region_img.ome.tif",
            cellprob_threshold = -1,
            flow_threshold     = 2,
            out_name           = "cp2_mask"
        )

