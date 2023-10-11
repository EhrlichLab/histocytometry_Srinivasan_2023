# install.packages("argparser")
suppressPackageStartupMessages(library(argparser, quietly= TRUE))

## Parse arguments 
parser <- arg_parser(description= "fcs")
parser <- add_argument(parser, "--path", type="character", help="path for csv segmentation quantification")
parser <- add_argument(parser, "--threshold", type="numeric", help="minimum cell size")
argv <- parse_args(parser)

## Filter area
df <- read.csv(argv$path)
print(dim(df))
df <- df[df$Area >= argv$threshold,]
print(dim(df))

write.csv(x= df, file= gsub(".csv", "_min_area_filtered.csv", argv$path), row.names= FALSE)