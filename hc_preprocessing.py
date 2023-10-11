## Run with Mesmer_pypi conda environment

import tifffile
import os 
import numpy as np
from skimage import ( data, restoration, util )
import pywt
import re
from scipy.ndimage import gaussian_filter
import typing
str_or_path = typing.Union[str, os.PathLike]

def preprocess_hc(hc_dir              : str_or_path,
                  image_dir           : str_or_path, 
                  rb_radius           : int, 
                  num_threads         : int, 
                  min_pixel_intensity : int,
                  gaussian_sigma      : int):

    combined_dir = os.path.join(hc_dir, image_dir)
    images       = os.listdir(combined_dir)
    image_name   = "image.tif"
    
    hc_image = tifffile.imread(os.path.join(combined_dir, image_name))
    num_slices = np.shape(hc_image)[0]
    filtered_images = [None] * num_slices
    
    if(rb_radius == 0):
        for slice in range(num_slices):
            slice_image= hc_image[slice, ...]

            ##-------------
            # Gaussian blur
            ##-------------
            blurred_slice = gaussian_filter(slice_image, 
                                            sigma= gaussian_sigma)
            filtered_images[slice] = blurred_slice

        filtered_image = np.stack(filtered_images)

        ##----------------
        # Writing gb image
        ##----------------
        out_dir = "".join(["gb", str(gaussian_sigma)])
        out_dir = os.path.join(hc_dir, "../processed_images", out_dir, image_dir)
    
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            ## os.makedirs needs to be used here instead of os.mkdir b/c I am making subdirectories. 
        tifffile.imwrite(os.path.join(out_dir, "image.ome.tif"), filtered_image)
        print("Done with an image")
        
    else:
        for slice in range(num_slices):
            slice_image= hc_image[slice, ...]   
        
            ##-------------
            # Gaussian blur
            ##-------------
            blurred_slice = gaussian_filter(slice_image,
                                            sigma= gaussian_sigma)
            ##----------------------
            # Background subtraction
            ##----------------------
            background = restoration.rolling_ball(image       = blurred_slice,
                                                  radius      = rb_radius,
                                                  num_threads = num_threads)
            ##------------------------
            # Writing background image
            ##------------------------
            back_out_name = "_".join(["background", str(slice), image_name])
            background_dir = os.path.join(combined_dir, "background")
            if not os.path.exists(background_dir):
                os.mkdir(background_dir)
            tifffile.imwrite(os.path.join(background_dir, back_out_name), background)

            filtered_images[slice] = blurred_slice - background 
            print("Done with a slice")

            ##---------------------------
            # Minimum threshold filtering
            ##---------------------------
            filtered_image = np.stack(filtered_images)
            filtered_image[filtered_image < min_pixel_intensity] = 0

        ##-------
        # Writing 
        ##-------
        out_dir = "".join(["gte", str(min_pixel_intensity), "_rb", str(rb_radius), "_gb", str(gaussian_sigma)])
            ## Saves how the images were processed 
        out_dir = os.path.join(hc_dir, "../processed_images", out_dir, image_dir)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            ## os.makedirs needs to be used here instead of os.mkdir b/c I am making subdirectories. 
        tifffile.imwrite(os.path.join(out_dir, "image.ome.tif"), filtered_image)
        print("Done with an image")


if __name__ == "__main__":
    raw_images_dir = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images"
    image_dirs = os.listdir(raw_images_dir)

    for image_dir in image_dirs:
        preprocess_hc(hc_dir              = raw_images_dir,
                      image_dir           = image_dir,
                      rb_radius           = 0, 
                      num_threads         = 70, 
                      min_pixel_intensity = 3,
                      gaussian_sigma      = 5)

