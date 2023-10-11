## Quick DC mask 
## Use with Mesmer_pypi conda environment
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
from sklearn.mixture import GaussianMixture
import tifffile
from pathlib import Path
import os
import cv2
import re
import typing
str_or_path = typing.Union[str, os.PathLike]
from deco import synchronized, concurrent


##---------
# Functions 
##---------

def setup_image(image_dir, glob_str= "image.tif", channel_loc= 2, median_ksize= 11, mean_ksize = (11,11)):
    image_name = [file_path for file_path in Path(image_dir).rglob(glob_str)]
    if len(image_name) > 1:
        print(image_dir + " has more than one tif")
    full_img  = tifffile.imread(os.path.join(image_dir, *image_name))
    img = full_img[channel_loc, ...]

    ## Normalize imaging data
    scaled_img = (img - np.mean(img)) / np.std(img)
    norm_img = np.arcsinh(scaled_img)

    return norm_img

def tissue_segmentation_plot(img_dir, img, thresholds, threshold_masks):
    print("Starting plot")
    _, ax = plt.subplots(1,len(thresholds)+1, figsize= (12,6))
    ax[0].imshow(img)
    ax[0].set_title("Marker")
    for i in range(1, len(thresholds)+1):
        ax[i].imshow(threshold_masks[i-1,...], cmap= plt.cm.gray)
        ax[i].set_title(f'Thresh: {round(thresholds[i-1], 1)}')
    plt.suptitle(img_dir) 
    plt.show() 

def run_image_GMM(img, n_gaussians= 4, downsample_divisor= 4):
    ## Downsampling to decrease computation time
    reshaped_img = img.reshape((img.size, 1)).flatten()

    ## Removing the outliers to deal w/ some bright medulla that might've caused some issues. 
        ## I thought about filtering the bottom outliers, but the values are always going to be 0. 
    top_outliers = np.percentile(reshaped_img, [97.5])
    reshaped_img = reshaped_img[(reshaped_img < top_outliers)]
    downsampled_img = np.random.choice(a       = reshaped_img, 
                                       size    = reshaped_img.size // downsample_divisor, 
                                       replace = False)

    ## Add dimension for use w/ GMM
    downsampled_img= downsampled_img[..., np.newaxis]
    gm = GaussianMixture(n_components= n_gaussians)
    gm.fit(downsampled_img)
        ## This is the step that takes a while 

    thresholds = gm.means_.flatten()
    thresholds.sort()

    return thresholds

def make_channel_threshold_masks(img, sorted_thresholds, remove_noise= False, min_spot_area= 100, median_ksize= 11):
    num_thresholds = len(sorted_thresholds)
    threshold_masks = np.zeros((num_thresholds, img.shape[0], img.shape[1]))

    for i in range(0, num_thresholds):
        threshold  = sorted_thresholds[i]
        binary_img = img > threshold
        binary_img = binary_img * np.uint8(1)     

        if remove_noise: 
            ## Small hole filling to capture missing nuclei
            contour, hier = cv2.findContours(binary_img, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
            for cnt in contour:
                cv2.drawContours(binary_img,[cnt],0,255,-1)

            ## Removing small noise
            binary_img = cv2.medianBlur(binary_img, median_ksize)
            num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(binary_img)

            ## Minimum area filtering
            filtered_img = np.zeros_like(binary_img)
            for label in range(1, num_labels):
                area = stats[label, cv2.CC_STAT_AREA]
                if area >= min_spot_area:
                    filtered_img[labels == label] = 255
        else:
            filtered_img = binary_img

        threshold_masks[i,...] = filtered_img
    return threshold_masks

@concurrent
def threshold_channel(img_dir, thresholds, channel_loc, out_name= "", glob_str= "reordered_image.ome.tif", save_mask= False, remove_noise= False, is_adaptive= False, show_plots= True) -> None:
    print("Starting " + img_dir + "\n")

    ## Take stain and apply thresholds
    img = setup_image(image_dir= img_dir, glob_str= glob_str,  channel_loc= channel_loc)
    
    if is_adaptive:
        threshold_masks = make_adaptive_binary_gate(img_arr = img) 
        show_plots = False
    else: 
        threshold_masks = make_channel_threshold_masks(img= img, sorted_thresholds= thresholds, remove_noise = remove_noise)
    threshold_masks = threshold_masks.astype("uint8") 

    if show_plots: 
        tissue_segmentation_plot(
            img_dir         = img_dir, 
            img             = img,
            thresholds      = thresholds, 
            threshold_masks = threshold_masks
        )

    if save_mask:
        assert min(threshold_masks.shape) == 1, "There are multiple threshold masks. Make sure there is only one threshold in the thresholds argument."
        assert len(out_name) > 0, "Set out_name to save image."
        tifffile.imwrite(os.path.join(img_dir, out_name), threshold_masks.squeeze())

    return None
    

@synchronized
def parallelize_channel_thresholding(img_dirs, thresholds, channel_loc, is_adaptive, save_mask, out_name= "", glob_str= "reordered_image.ome.tif", remove_noise= False, show_plots= True):
    ## Use this function to speed up plotting of whole lobe masks and medulla masks
    ## I'll visually inspect the results to pick the best thresholds
    for img_dir in img_dirs:
        threshold_channel(
            img_dir      = img_dir, 
            thresholds   = thresholds, 
            channel_loc  = channel_loc, 
            glob_str     = glob_str,
            save_mask    = save_mask,
            out_name     = out_name, 
            remove_noise = remove_noise,
            is_adaptive  = is_adaptive,
            show_plots   = show_plots
        )

def overlay_mask(img_dir, image_name, mask_name, out_name) -> None:
    img= tifffile.imread(os.path.join(img_dir, image_name))
    mask = tifffile.imread(os.path.join(img_dir, mask_name))
    mask = mask[np.newaxis,...]

    img_w_mask= np.concatenate((img, mask))

    tifffile.imwrite(os.path.join(img_dir, out_name), img_w_mask)

    return None

def make_adaptive_binary_gate(img_arr, median_kernel_size= 5, adaptive_subtraction= 1, adaptive_kernel_size= 5, min_spot_area= 0):
    ## OpenCV documentation recommended blurring prior to adaptive thresholding
    img_arr = img_arr.astype("uint8")
    blurred = img_arr
    blurred = cv2.medianBlur(img_arr, median_kernel_size)

    threshold_img = cv2.adaptiveThreshold(blurred, 
                                          255, 
                                          cv2.ADAPTIVE_THRESH_GAUSSIAN_C, ## I think Gaussian looks significantly better than mean for cortical images.
                                          cv2.THRESH_BINARY, 
                                          adaptive_kernel_size, 
                                          adaptive_subtraction)
    filtered_img = threshold_img[np.newaxis, ...]



    # ## Small hole filling to capture missing nuclei
    # contour, hier = cv2.findContours(threshold_img, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    # for cnt in contour:
    #     cv2.drawContours(threshold_img,[cnt],0,255,-1)

    # ## Filtering to remove salt and pepper noise
    # threshold_img = cv2.medianBlur(threshold_img, median_kernel_size)

    # ## Removing specks
    # num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(threshold_img)
    # filtered_img = np.zeros_like(threshold_img)
    # for label in range(1, num_labels):
    #     area = stats[label, cv2.CC_STAT_AREA]
    #     if area >= min_spot_area:
    #         filtered_img[labels == label] = 255

    # filtered_img = filtered_img[np.newaxis, ...]
    
    return filtered_img


def segment_image_using_mask(img_dir, image_name= "reordered_image.ome.tif", mask_name= "DC_mask_adaptive.tif", out_name= "DC_image.ome.tif") -> None:

    img = tifffile.imread(os.path.join(img_dir, image_name))
    DC_binary_mask = tifffile.imread(os.path.join(img_dir, mask_name))

    DC_sub_img = img * DC_binary_mask

    tifffile.imwrite(os.path.join(img_dir, out_name), DC_sub_img)

    return None 



if __name__=="__main__":

    raw_dir = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images/images_2023-08-10"
    img_dirs = [os.path.join(raw_dir, img_dir) for img_dir in os.listdir(raw_dir) if re.search("_[A-D]$", img_dir)]

    parallelize_channel_thresholding(
        img_dirs     = img_dirs,
        thresholds   = np.arange(1.75, 2, 0.25),
        channel_loc  = 2,
        save_mask    = True,
        out_name     = "DC_mask.tif",
        is_adaptive  = False,
        remove_noise = True,
        show_plots   = False
    )
    ## 1.5 was too large, 2 was too small, maybe 1.75 will be the answer. 
    ## I need to look at this mask to see how it looks. If it's at least close, I'll take it. 
    ## I lowered the threshold and increased the area filtering. 
