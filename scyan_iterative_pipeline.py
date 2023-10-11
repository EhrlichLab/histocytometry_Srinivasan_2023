##---------------
# Scripting scyan
##---------------
import scyan
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import anndata
import re 
from pathlib import Path
## Run w/ scyan conda environment 



##---------
# Functions
##---------
def write_knowledge_tables(out_dir):
     ## This function will write the hardcoded knowledge_tables and will serve as a reference. 
        ## I can also quickly iterate and re-write the tables. 
     DAPI_CD63_CD11c_Sirpa_table = pd.DataFrame(data= {"CD63"        : [1, -1],
                                                       "CD11c"       : [np.nan, 1], ## CD11c turned down somewhat for both (I don't know if Scyan will care about this.)
                                                       "Sirpa"       : [1, 1]}, ## aDC2 turns down Sirpa 
                                                index= ["aDC2", "cDC2"])
     DAPI_CD63_CD11c_Sirpa_table.to_csv(os.path.join(out_dir, "DAPI_CD63_CD11c_Sirpa_table.csv"), index= True)


     DAPI_CD63_CD11c_XCR1_table = pd.DataFrame(data= {"CD63"  : [1, -1],
                                                      "CD11c" : [np.nan, 1],
                                                      "XCR1"  : [1, 1]}, ## XCR1 downregulated on aDC1. 
                                               index= ["aDC1", "aDC2"])
     DAPI_CD63_CD11c_XCR1_table.to_csv(os.path.join(out_dir, "DAPI_CD63_CD11c_XCR1_table.csv"), index= True)


     DAPI_Sirpa_CD11c_MerTK_table = pd.DataFrame(data= {"Sirpa" : [1],
                                                        "CD11c" : [np.nan],
                                                        "MerTK" : [1]},
                                                index= ["macrophage"])
     DAPI_Sirpa_CD11c_MerTK_table.to_csv(os.path.join(out_dir, "DAPI_Sirpa_CD11c_MerTK_table.csv"), index= True)
    

     DAPI_Sirpa_CD11c_CD14_table = pd.DataFrame(data= {"Sirpa" : [1],
                                                       "CD11c"  : [1],
                                                       "CD14"   : [1]},
                                                index= ["cDC2"])
     DAPI_Sirpa_CD11c_CD14_table.to_csv(os.path.join(out_dir, "DAPI_Sirpa_CD11c_CD14_table.csv"), index= True)

     DAPI_B220_CD11c_SiglecH_table = pd.DataFrame(data= {"B220"    : [1],
                                                         "CD11c"   : [np.nan], #intermediate
                                                         "SiglecH" : [1]},
                                                  index= ["pDC"])
     DAPI_B220_CD11c_SiglecH_table.to_csv(os.path.join(out_dir, "DAPI_B220_CD11c_SiglecH_table.csv"), index= True)

     DAPI_CD207_CD11c_XCR1_table = pd.DataFrame(data= {"CD207"    : [1],
                                                       "CD11c"   : [np.nan], #intermediate
                                                       "XCR1"     : [1]},
                                                  index= ["cDC1"])
     DAPI_CD207_CD11c_XCR1_table.to_csv(os.path.join(out_dir, "DAPI_CD207_CD11c_XCR1_table.csv"), index= True)


     ## Write dictionary output 
     img_experiments = {"DAPI_CD63_CD11c_Sirpa"   : os.path.join(out_dir, "DAPI_CD63_CD11c_Sirpa_table.csv"), 
                        "DAPI_CD63_CD11c_XCR1"    : os.path.join(out_dir, "DAPI_CD63_CD11c_XCR1_table.csv"), 
                        "DAPI_B220_CD11c_SiglecH" : os.path.join(out_dir, "DAPI_B220_CD11c_SiglecH_table.csv"), 
                        "DAPI_Sirpa_CD11c_CD14"   : os.path.join(out_dir, "DAPI_Sirpa_CD11c_CD14_table.csv"), 
                        "DAPI_Sirpa_CD11c_MerTK"  : os.path.join(out_dir, "DAPI_Sirpa_CD11c_MerTK_table.csv")}
     return(img_experiments)

## I might include CD11c-positivity as a filtering step. I'm not sure if I want to gate on CD11c expression prior to Scyan and then call cells based on the other variables or if I want to call all cells as CD11c+ or NA. 


def combine_image_data(image_dir, img_exp_type):
    csv_paths      = [csv_path for csv_path in Path(image_dir).rglob('*_mask.csv*') if re.search(img_exp_type, str(csv_path))]
        ## i.e. take all of the cells from the DAPI_CD63_CD11c_Sirpa images and combine them for cell type identification. 
        ## This assumes no batch differences between images, but it means that the same cell type is identified in the same way across all samples. 
    mask_df_list   = [None] * len(csv_paths)

    for i in range(0, len(csv_paths)):
        csv_path = csv_paths[i]
        mask_df_list[i] = pd.read_csv(csv_path)
        mask_df_list[i]["experiment"] = str(csv_path)

    combined_mask_df = pd.concat(mask_df_list, axis= 0)
    codes, uniques = pd.factorize(combined_mask_df.experiment)
        ## uniques is what I'll need to identify which cluster is which experiment later
        ## I have to do this because scyan can only handle numeric metadata values. 
    combined_mask_df.experiment = codes
    
    data_dir= os.path.join(image_dir, "data")
    if not os.path.exists(data_dir):
            os.makedirs(data_dir)
    combined_mask_path = os.path.join(data_dir, "combined_mask_df_" + img_exp_type + ".csv")
    combined_mask_df.to_csv(combined_mask_path, index= False)
    experiment_dummy = pd.DataFrame({"experiment" : np.unique(codes),
                                     "exp_name"    : uniques})
    experiment_dummy.to_csv(os.path.join(data_dir, "experiment_dummy_" + img_exp_type + ".csv"), index= False)

    return(combined_mask_df)


def convert_df_to_adata(df, marker_list):
    adata = anndata.AnnData(df)
    adata.obs = adata[:, ["CellID", "X_centroid", "Y_centroid", "experiment"]].to_df()
    adata = adata[:, marker_list]
                       
    print(f"Created anndata object with {adata.n_obs} cells and {adata.n_vars} markers.\n\n-> The markers names are: {', '.join(adata.var_names)}\n-> The non-marker names are: {', '.join(adata.obs.columns)}")

    return(adata)

def preprocess_scyan(adata, is_cytof):
    if is_cytof:
        scyan.preprocess.asinh_transform(adata)
    else: #
        scyan.preprocess.auto_logicle_transform(adata)
            ## Some transformation designed for flow: https://pubmed.ncbi.nlm.nih.gov/16604519/
    scyan.preprocess.scale(adata)
    return(adata)


def full_scyan_pipeline(image_dir, 
                        knowledge_table_path, 
                        img_exp_type, 
                        num_workers = 1, 
                        is_cytof    = False):
     ## num_workers doesn't work if their multiple processes for some reason. 

     ## Making data
     knowledge_table = pd.read_csv(knowledge_table_path)

     combined_mask_df = combine_image_data(image_dir = image_dir, img_exp_type = img_exp_type)
     adata = convert_df_to_adata(df= combined_mask_df, marker_list= knowledge_table.columns.to_list())
     adata = preprocess_scyan(adata = adata, is_cytof = is_cytof)


     ## Running Scyan
     model = scyan.Scyan(adata, knowledge_table)
     model.fit(num_workers= num_workers)
     model.predict()


     ## Plotting data
     scyan.tools.umap(adata, markers= knowledge_table.columns)
     scyan.plot.umap(adata, color = "scyan_pop", 
                     title        = "Scyan predictions")
     plt.savefig(os.path.join(image_dir, "plots/scyan_plots", img_exp_type + "_scyan_umap_plot.jpg"))
     
     print(adata.obs.scyan_pop.value_counts())
     print(adata.obs.scyan_pop.value_counts(normalize= True))


     scyan.plot.scatter(adata, 
                        population = None, 
                        markers    = knowledge_table.columns.to_list(),
                        show       = False) ## Could this be colored by population?
     plt.savefig(os.path.join(image_dir,"plots/scyan_plots", img_exp_type + "_scyan_scatter_plot.jpg"))   


     scyan.plot.pops_expressions(model, 
                                 latent  = True, 
                                 figsize = (10, 6), 
                                 show    = False)
     plt.savefig(os.path.join(image_dir,"plots/scyan_plots", img_exp_type + "_scyan_pop_expression.jpg"))   


     ## Write output
     adata.obs.to_csv(os.path.join(image_dir, "data", img_exp_type + "_combined_mask_df_scyan_obs.csv"), index= False)
     pd.DataFrame(adata.X, columns = adata.var_names).to_csv(os.path.join(image_dir, "data", img_exp_type + "_combined_mask_df_scyan_X.csv"),index= False)



if __name__ == "__main__":
    ## For each experiment set-up, we expect to find different cell types and we have different markers. 
    ## Scyan pulls this data from the knowledge table, so I just need to subset the img_dirs and select the knowledge_table that is relevant for the experiment that
    ## I am analyzing. 

    ## Right now, this is working with the max mask. Doing this with a different mask would require either having a different top direcotory or naming the .csv path 
    ## in a way that holds the max mask data. 


    img_path = "/stor/scratch/Ehrlich/Users/John/histocytometry/raw_images/images_2023-08-10"
    img_experiments = write_knowledge_tables(out_dir= os.path.join(img_path, "knowledge_tables"))

    for img_experiment in img_experiments.keys(): 
        full_scyan_pipeline(image_dir            = img_path, 
                            knowledge_table_path = img_experiments[img_experiment],
                            img_exp_type         = img_experiment)
        ## I'm not super pleased with how the pDC gating looks. 
        ## It could be that Scyan thinks the cells it's omitting are outliers or part of another latent distribution. 
        ## It could also be a normalization issue. 
      
