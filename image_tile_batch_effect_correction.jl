## I want to do the Lau Lab's tile normalization just to see if it helps a little bit. 

## I'll do mean divsion followed by log10 transformation to make the data less skewed. 
using FileIO
using Images
using StatsBase
# using OMETIFF


##---------
# Functions 
##---------
function slide_batch_effect_mean_divide_log10(img :: Array{T, 3}, tile_length :: Int) where T 
    ## This type hinting says that img must be a three dimensional array of any type (so long as the type is consistent)

    ## Assumes square tile size  
    img_size = size(img)

    scaled_image = zeros(img_size)
    num_x_tiles = img_size[1] ÷ tile_length
    num_y_tiles = img_size[2] ÷ tile_length
        ## Julia loves unicode representations for operators.
        ## This is the integer divide. 

    for channel_num in 1:img_size[3]        
        for i in 1:num_x_tiles 
            start_x = (i-1) * tile_length + 1
            end_x   = i * tile_length
            
            for j in 1:num_y_tiles 
                start_y = (j-1) * tile_length + 1
                end_y   = j * tile_length
            
                img_tile = img[start_x:end_x, start_y:end_y, channel_num]
            
                channel_mean = mean(img_tile)

                scaled_image[start_x:end_x, start_y:end_y, channel_num] = log10.(img_tile ./ channel_mean)
            end
        end
    end

    
    return scaled_image
end

function batch_correction_wrapper(img_file :: AbstractString, out_file :: AbstractString)
    @assert isfile(img_file) "Image file doesn't exist or isn't a file."

    ## Load image & convert to original data type 
    img = rawview(channelview(load(img_path)))
    new_img = slide_batch_effect_mean_divide_log10(img, tile_size)

    ## Save image (have to convert back to image datatype)
    scaled_image = reinterpret(Gray{eltype(scaled_image)}, scaled_image)
    save(img, out_file) 
end


##----------------
# Preparing images
##---------------- 
raw_dir = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images"
img_dirs = [joinpath(raw_dir, img_dir) for img_dir in readdir(raw_dir) if contains(img_dir, "20x") & contains(img_dir, r"_(A|B|C)$")]

imgs_to_replace= joinpath.(raw_dir, ["20x_pan_DAPI_B220_CD11c_SiglecH_A", 
                                     "20x_pan_DAPI_B220_CD11c_SiglecH_B", 
                                     "20x_pan_DAPI_CD207_CD11c_XCR1_B",
                                     "20x_pan_DAPI_CD207_CD11c_XCR1_C",
                                     "20x_pan_DAPI_Sirpa_CD63_MerTK_C"])

img_dir = [img_dir for img_dir in img_dirs if !(img_dir in imgs_to_replace)]

tile_size= 460
    ## This is a placeholder. 


##------------------
# Normalizing images 
##------------------
for img_dir in img_dirs
    img_file = joinpath(img_dir, "image.tif")
    out_file = joinpath(img_dir,  "mean_log10_image.tif")

    batch_correction_wrapper(img_file, out_file)                               
end


##---------------------------------------
# Validating batch effect correction code
##---------------------------------------
base_array = ones(3,3)
test_channel = vcat(hcat(base_array .* 1, base_array .* 2, base_array .* 3),
                    hcat(base_array .* 4, base_array .* 5, base_array .* 6),
                    hcat(base_array .* 7, base_array .* 8, base_array .* 9))
test_image = cat(test_channel, test_channel, test_channel, dims= 3)

actual_output = slide_batch_effect_mean_divide_log10(test_image, 3)
expected_output= zeros(9,9,3)

all(expected_output .== actual_output)
    ## Returns true

## THe starting tile doesn't have the overlap, so maybe I need to have an initial tile setting for histocytometry. 
## I need the percent overlap 
## Lauren worries that the bright spot may be due to a lack of changing filter on the microscope. 


