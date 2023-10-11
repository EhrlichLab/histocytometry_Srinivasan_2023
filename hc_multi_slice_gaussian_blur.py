## Run with Mesmer_pypi conda environment

import tifffile
import os 
import numpy as np
from skimage import ( data, restoration, util )
import pywt
import re
from scipy.ndimage import gaussian_filter

hc_dir     = "../"
image_dirs    = [image_dir for image_dir in os.listdir(hc_dir) if not re.match("\.", image_dir) and re.match("40x", image_dir)]

image_dirs= ["noBack_cortex_radius","not a dir"]
for image_dir in image_dirs:
    images     = os.listdir(os.path.join(hc_dir, image_dir))
    image_name =  "".join([image for image in images if re.search(".tif", image) and not re.search("ack|mask|gating|.ome.tif|blur", image)])

    hc_image = tifffile.imread(os.path.join(hc_dir, image_dir, image_name))
    num_rows = np.shape(hc_image)[0]
    filtered_images = [None] * num_rows

    for slice in range(num_rows):
        slice_image = hc_image[slice, ...]
        new_slice   = gaussian_filter(slice_image, sigma= 5)
        filtered_images[slice] = new_slice
        print("Done with a slice")
    
    
    filtered_image = np.stack(filtered_images)
    out_name = "_".join(["blur", image_name])
    out_name = re.sub(".tif|.tiff", ".ome.tif", out_name)
    tifffile.imwrite(os.path.join(hc_dir, image_dir, out_name), filtered_image)

    print("Done with an image")
