## Run with Mesmer_pypi conda environment

import tifffile
import os 
import numpy as np
from skimage import ( data, restoration, util )
import pywt
import re

hc_dir     = "../"
# image_dirs    = [image_dir for image_dir in os.listdir(hc_dir) if re.match("40x", image_dir) and re.search("Sirpa", image_dir) and not re.search("pan", image_dir)]
image_dirs = ["noBack_cortex_radius", "not_a_dir"] 

for image_dir in image_dirs:
    images     = os.listdir(os.path.join(hc_dir, image_dir))
    image_name =  "".join([image for image in images if re.search(".tif", image) and not re.search("sub|mask|.ome.tif", image)])

    hc_image = tifffile.imread(os.path.join(hc_dir, image_dir, image_name))
    num_rows = np.shape(hc_image)[0]
    filtered_images = [None] * num_rows

    for slice in range(num_rows):
        slice_image= hc_image[slice, ...]
        background = restoration.rolling_ball(image       = slice_image,  
                                              radius      = 100,
                                              num_threads = 70)
        back_out_name = "_".join(["background", str(slice), image_name])
        ## Save background subtraction reference to another subdirectory
        ## For minimum threshold, would I want to get rid of anything below an 
        ## absolute pixel intensity, or should I set a percentile threshold?
        ## If I do a percentile threshold, it should be of non-zero pixels.
        ## The issue of an absolute threshold is that not all images have the same intensity. 
        tifffile.imwrite(os.path.join(hc_dir, image_dir, back_out_name), background)
        filtered_images[slice] = slice_image - background 

        print("Done with a slice")
    
    
    filtered_image = np.stack(filtered_images)
    out_name = "_".join(["noBack", image_name])
    out_name = re.sub(".tiff|.tif", ".ome.tif", out_name)
    tifffile.imwrite(os.path.join(hc_dir, image_dir, out_name), filtered_image)

    print("Done with an image")
