import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import tifffile
from pathlib import Path
import os
import cv2
import re
import typing
str_or_path = typing.Union[str, os.PathLike]
from deco import concurrent, synchronized


@concurrent
def GMM_wrapper(arr, n ):
    ## Making this wrapper to add deco concurrency to GMM function
    model = GaussianMixture(n_components= n).fit(arr)
    return model

@synchronized # And we add this for the function which calls the concurrent function
def parallelize_GMMs(arr, max_components):
  models = [GMM_wrapper(arr= arr, n= n) for n in np.arange(2, max_components)]
  return(models)




##--------
# Test run
##--------


## Parameters
image_dir= "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images/20x_pan_DAPI_CD207_CD11c_XCR1_A"
median_ksize       = 11 ## has to be odd
morph_ksize        = (75,75) ## as opposed to (200,200)
final_median_ksize = 251 # Bigger to smooth out medulla regions
channel_loc        = 2
out_name           = "medulla"
inv_name           = "cortex_and_background"
glob_str           = 'image.tif'


## Loading image
image_name = [file_path for file_path in Path(image_dir).rglob(glob_str)]
if len(image_name) > 1:
    print(image_dir + " has more than one tif")
full_img  = tifffile.imread(os.path.join(image_dir, *image_name))
img = full_img[channel_loc, ...]

## Median blur
med_blur = cv2.medianBlur(img, median_ksize)

## Trying a mean blur instead to get an even blurry medulla.
mean_ksize = (11,11)
mean_blur = cv2.blur(img, mean_ksize)

downsampled_img = np.random.choice(a       = img.reshape((mean_blur.size, 1)).flatten(),
                                   size    = img.size // 4,
                                   replace = False)

## Add dimension for use with GMM
downsampled_img= downsampled_img[..., np.newaxis]


n_components = np.arange(2,10)
models= parallelize_GMMs(arr= downsampled_img, max_components= 10)

import pickle 

## wb means "write binary"
## listfile is the saved name 
## fp is the open "file pickle" 
pickle_file= os.path.join(image_dir, "model_list.pkl")
with open(pickle_file, 'wb') as fp:
    pickle.dump(models, fp)
## The dump is for dumping the models list into the oepn file.
## The with indentation sets the boundaries for where the file is open.

## I'll take this model list into my .ipynb for the image comparisons. 


