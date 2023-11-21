# title: "Joining gated data to new CMJ region"
# author: John Moore
# date: 2023-09-18

library(dplyr)
library(data.table)
library(here)
source("/stor/work/Ehrlich/Users/John/projects/AKOYA/scripts/AKOYA_scripts/scripts/AKOYA_clustering_comparison_functions.R")


##---------
# Functions
##---------
fcs_to_df= function(fcs_path){
  df <- as.data.frame(exprs(read.FCS(fcs_path, truncate_max_range= FALSE)))
  return(df)
}

combine_fcs_output_files <- function(histo_dir, pos_keyword, neg_keyword, mask_quant_keyword){
  pop_files = list.files(path       = histo_dir, 
                         pattern    = pos_keyword,
                         full.names = TRUE)
  pop_files = pop_files[!grepl(neg_keyword, pop_files)]
  print(paste("Files used", pop_files))
  
  len_files = length(pop_files)
  dfs_list = vector(mode= "list", length= len_files)
  
  ## Re-format flow data w/ relevant metadata
  for(i in 1:len_files){
    hc_file = pop_files[i]
    tmp_df = fcs_to_df(hc_file)
    tmp_df$cell_type = gsub(paste0(histo_dir, "/|export_|null_|.fcs"), "", hc_file)
    dfs_list[[i]] = tmp_df
  }
  combo_df = do.call(rbind, dfs_list) %>% 
    select(CellID, cell_type) 
  
  ## Pull in predicted cells and add cell type annotation
  mask_file = list.files(path = histo_dir, 
                         pattern = mask_quant_keyword,
                         full.names = TRUE)
  mask_df = as.data.frame(fread(mask_file))
  out_df= left_join(mask_df, combo_df, by= "CellID") %>% 
    mutate(cell_type = ifelse(is.na(cell_type), "not_annotated", cell_type)) 
  
  return(out_df)
}




##--------------
# Running images
##--------------
raw_dir = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images/images_2023-08-10"
hc_dirs = list.dirs(path= raw_dir, full.names= TRUE) 
hc_dirs = sort(hc_dirs[grepl("_[A-D]$", hc_dirs)])
hc_dirs = hc_dirs[!grepl("Sirpa_C$|CD207_CD11c_XCR1_C$|MerTK_C$", hc_dirs)]



for(histo_dir in hc_dirs){
  print(gsub(raw_dir, "", histo_dir))
  all_gated_pops <- combine_fcs_output_files(histo_dir= histo_dir, 
                                             pos_keyword= "export_null_.*.fcs$",
                                             neg_keyword= "min_area_filtered|.wsp",
                                             mask_quant_keyword = "reordered_image_w_thymus_regions_50um_cp2_mask_min_area_filtered.csv") %>% 
    select(-contains("tissue"))

  tissue_regions <- paste0(histo_dir, "/reordered_image_w_thymus_regions_50um_cp2_mask_min_area_filtered.csv") %>% 
    fread() %>% 
    as.data.frame() %>% 
    select(CellID, contains("tissue"))
  
  print(all(all_gated_pops$CellID %in% tissue_regions$CellID))
  gated_pops_w_new_CMJ <- left_join(all_gated_pops, tissue_regions, by= c("CellID"))
  
  ## I need to split ties in the median tissue value.
    ## I'm going to split the tie randomly.
  ties_df <- gated_pops_w_new_CMJ[gated_pops_w_new_CMJ$tissue_intensity_median %% 1 != 0,]
  ties_df <- ties_df %>% 
    mutate(rank = sample(c(0,1),size = nrow(.), replace= TRUE),
           tissue_intensity_median = ifelse(rank == 0, 
                                            ceiling(tissue_intensity_median), 
                                            floor(tissue_intensity_median))) %>% 
    select(-rank)
  gated_pops_w_new_CMJ <- gated_pops_w_new_CMJ[gated_pops_w_new_CMJ$tissue_intensity_median %% 1 == 0,]
  gated_pops_w_new_CMJ <- rbind(ties_df, gated_pops_w_new_CMJ)
    
    
  print(dim(gated_pops_w_new_CMJ))
  fwrite(gated_pops_w_new_CMJ,
         paste0(histo_dir, "/", "gated_pops_w_new_CMJ.csv"),
         row.names= FALSE)
}


