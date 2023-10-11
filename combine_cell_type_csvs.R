options(install.packages.compile.from.source = "always")
packages= c("argparser")
install.packages(setdiff(packages, rownames(installed.packages())), type= "both")
  ## Installing any packages that aren't found

## Parse arguments 
parser <- arg_parser(description= "fcs")
parser <- add_argument(parser, "--path", type="character", help="path for csv segmentation quantification")
argv   <- parse_args(parser)


## Combine cell type csvs
cell_type_files = list.files(path= argv$path, full.names= TRUE, recursive= TRUE, pattern= "_mask_area_filtered.csv")

cell_type_dfs= vector(mode= "list", length= length(cell_type_files))
for(i in 1:length(cell_type_files)){
    cell_type_dfs[[i]]= read.csv(cell_type_files[i])
}

combined_cell_type_df = do.call(rbind, cell_type_dfs)

write.csv(x= combined_cell_type_df, x= paste0(argv$path,"/combined_cell_type_masks_df.csv"), row.names= FALSE)


