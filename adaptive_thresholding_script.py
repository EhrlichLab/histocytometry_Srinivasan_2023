## I am re-writing adaptive_thresholding, so I can script it in the background. 
## Use Mesmer_pypi conda environment
import numpy as np
import cv2
import os
import matplotlib.pyplot as plt
import tifffile
import pandas as pd
import re


##----------
# Functions
##----------
def make_adaptive_binary_gate(img_arr, median_kernel_size= 21, adaptive_subtraction= 0, adaptive_kernel_size= 21, morphology_kernel = (15,15)):
    num_channels = np.shape(img_arr)[0] 
    ## tifffile loads in ch x h x w

    threshold_list = [None] * (num_channels -1)
        ## Removing DAPI from morphological masks b/c not used for this segmentation

    for i in range(1, num_channels):
        slice = img_arr[i,...]
        
        ## OpenCV documentation recommended blurring prior to adaptive thresholding
        blurred = cv2.medianBlur(slice, median_kernel_size)

        threshold_img = cv2.adaptiveThreshold(blurred, 
                                              255, 
                                              cv2.ADAPTIVE_THRESH_GAUSSIAN_C, ## I think Gaussian looks significantly better than mean for cortical images.
                                              cv2.THRESH_BINARY, 
                                              adaptive_kernel_size, 
                                              adaptive_subtraction)
        
        ## Morphological operation to remove noise
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, morphology_kernel)
        morphology_img = cv2.morphologyEx(threshold_img, cv2.MORPH_CLOSE, kernel)

        ## Adding channel dimension for concatenation 
        threshold_list[i-1] = morphology_img[np.newaxis, ...]

    adaptive_binary_img = np.concatenate(threshold_list)
    return(adaptive_binary_img)


def make_cell_type_mask(adaptive_binary_img, cell_type):
    img_shape = np.shape(adaptive_binary_img)[1:3]
    mask = np.zeros(img_shape)

    ## Channel gates assume that channel order matches directory names in /stor/scratch/Ehrlich/Users/John/histocytometry/raw_images 
        ## See also /stor/work/Ehrlich/Users/John/projects/misc/histocytometry/scripts/scyan_iterative_pipeline.py 
        ## No DAPI channel for gates 
    if cell_type == "aDC2_Sirpa":
        for row in range(img_shape[0]):
            for col in range(img_shape[1]):
                pixel_val = adaptive_binary_img[:, row, col]
                mask[row, col] = np.all(pixel_val == [255,0,255]) | np.all(pixel_val == [255,255,255])
        print(pixel_val)
    elif cell_type == "cDC2":
        for row in range(img_shape[0]):
            for col in range(img_shape[1]):
                pixel_val = adaptive_binary_img[:, row, col]
                mask[row, col] = np.all(pixel_val == [0,255,255])
        print(pixel_val)
    elif cell_type == "aDC1":
        for row in range(img_shape[0]):
            for col in range(img_shape[1]):
                pixel_val = adaptive_binary_img[:, row, col]
                mask[row, col] = np.all(pixel_val == [255,0,255]) | np.all(pixel_val == [255,255,255])
    elif cell_type == "aDC2_XCR1":
        for row in range(img_shape[0]):
            for col in range(img_shape[1]):
                pixel_val = adaptive_binary_img[:, row, col]
                mask[row, col] = np.all(pixel_val == [0,255,255])
    elif cell_type == "macrophage":
        for row in range(img_shape[0]):
            for col in range(img_shape[1]):
                pixel_val = adaptive_binary_img[:, row, col]
                mask[row, col] = np.all(pixel_val == [255,0,255])
    elif cell_type == "monoDC2":
        for row in range(img_shape[0]):
            for col in range(img_shape[1]):
                pixel_val = adaptive_binary_img[:, row, col]
                mask[row, col] = np.all(pixel_val == [255,255,255])
    elif cell_type == "pDC":
        for row in range(img_shape[0]):
            for col in range(img_shape[1]):
                pixel_val = adaptive_binary_img[:, row, col]
                mask[row, col] = np.all(pixel_val == [255, 0, 255]) | np.all(pixel_val == [255,255,255])
    else: 
        raise ValueError("Cell type: " + cell_type + " doesn't match hardcoded options")
    
    return(mask)    

def np_table(arr):
    ## Making this to get some summary stats on the images conveniently
    unique_values, counts = np.unique(arr, return_counts=True)
    np_table = pd.DataFrame({"unique_values" : unique_values,
                            "counts" : counts})
    return(np_table)
    

def gate_thymus_histocytometry(img_dir, img_name, cell_type, median_kernel_size= 21):
    ## I should think about doing kwargs for make_adaptive_binary_img parameters.

    ## Making individual marker binary gates and saving them for later reference. 
    img_arr = tifffile.imread(os.path.join(img_dir, img_name))
    adaptive_binary_img = make_adaptive_binary_gate(img_arr, median_kernel_size= median_kernel_size)
    out_path = os.path.join(img_dir, "individual_adaptive_threshold_gates.ome.tif")
    tifffile.imwrite(out_path, adaptive_binary_img)
    print("Done with adaptive thresholding!")

    ## Making cell type mask using combined binary mask input 
    cell_type_mask = make_cell_type_mask(adaptive_binary_img, cell_type)
    cell_type_mask = cell_type_mask * 255
        ## Converting 0 and 1 to more dynamic color range.
    cell_type_mask = cell_type_mask.astype("uint8")
    print("Done with cell type mask!")

    print(np_table(cell_type_mask))

    ## Make cell type directory for easier use with MCQuant later.
    out_path = os.path.join(img_dir, "cell_type_" + cell_type)
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    tifffile.imwrite(os.path.join(out_path, cell_type + "_mask.tif"), cell_type_mask)

        ## I'll run this through MCquant to get the cell count and area 
        ## I can use the area to filter out garbage "cells" and count easier than doing it by hand. 
            ## There may be an easier way to do this with FIJI.



if __name__ == '__main__':
    
    ## Preparing data
    raw_dir  = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images"
    img_dirs = [os.path.join(raw_dir, my_dir) for my_dir in os.listdir(raw_dir)]

    DAPI_CD63_CD11c_Sirpa_dirs    = [img_dir for img_dir in img_dirs if re.search("DAPI_CD63_CD11c_Sirpa",   img_dir)]
    DAPI_CD63_CD11c_XCR1_dirs     = [img_dir for img_dir in img_dirs if re.search("DAPI_CD63_CD11c_XCR1",    img_dir)]
    DAPI_CD63_CD11c_MerTK_dirs    = [img_dir for img_dir in img_dirs if re.search("DAPI_CD63_CD11c_MerTK",   img_dir)]
        ## I don't actually have any CD63_CD11c_MerTK images.
    DAPI_Sirpa_CD11c_CD14_dirs    = [img_dir for img_dir in img_dirs if re.search("DAPI_Sirpa_CD11c_CD14",   img_dir)]
    DAPI_B220_CD11c_SiglecH_dirs  = [img_dir for img_dir in img_dirs if re.search("DAPI_B220_CD11c_SiglecH", img_dir)]

    ## Temporarily removing dirs that I've already run
    if False: 
        ## Running images
        for indiv_dir in DAPI_CD63_CD11c_Sirpa_dirs:
            if re.search("pan", indiv_dir):
                print(indiv_dir)
                gate_thymus_histocytometry(img_dir = indiv_dir, img_name = "sub_image.ome.tif", cell_type = "aDC2_Sirpa", median_kernel_size= 31)
                gate_thymus_histocytometry(img_dir = indiv_dir, img_name = "sub_image.ome.tif", cell_type = "cDC2",       median_kernel_size= 31)

        for indiv_dir in DAPI_CD63_CD11c_XCR1_dirs:
            if re.search("pan", indiv_dir):
                print(indiv_dir)
                gate_thymus_histocytometry(img_dir = indiv_dir, img_name = "sub_image.ome.tif", cell_type = "aDC1")
                gate_thymus_histocytometry(img_dir = indiv_dir, img_name = "sub_image.ome.tif", cell_type = "aDC2_XCR1")
            ## Is this also aDC2? 

        for indiv_dir in DAPI_CD63_CD11c_MerTK_dirs:
            if re.search("pan", indiv_dir):
                print(indiv_dir)
                gate_thymus_histocytometry(img_dir = indiv_dir, img_name = "sub_image.ome.tif", cell_type = "macrophage")

    for indiv_dir in DAPI_Sirpa_CD11c_CD14_dirs:
        if re.search("pan", indiv_dir):
            print(indiv_dir)
            gate_thymus_histocytometry(img_dir = indiv_dir, img_name = "sub_image.ome.tif", cell_type = "monoDC2")

    for indiv_dir in DAPI_B220_CD11c_SiglecH_dirs:
        if re.search("pan", indiv_dir):
            print(indiv_dir)
            gate_thymus_histocytometry(img_dir = indiv_dir, img_name = "sub_image.ome.tif", cell_type = "pDC")



## Loop through raw dirs
## Loop through cell type directories for each image. 
## Run MCQuant on each image and cell type to extract the data 
## Combine cell type csv files 
## Filter out cells that are below a certain size 
## Overlay the masks onto the original image. 
