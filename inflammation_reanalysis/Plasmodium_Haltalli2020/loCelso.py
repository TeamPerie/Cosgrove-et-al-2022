import os
import scanpy as sc
import anndata as ad
import csv
import pandas
adata.obs


os.chdir("/Users/jasoncosgrove/Desktop/wetransfer-2a31f3/")

adata = ad.read_h5ad("malaria_control_data_genespector.h5ad")
adata.var_names_make_unique()
adata.obs.to_csv('control_metadata.csv')

adata = ad.read_h5ad("malaria_infected_data_genespector.h5ad")
adata.var_names_make_unique()
adata.obs.to_csv('infected_metadata.csv')
