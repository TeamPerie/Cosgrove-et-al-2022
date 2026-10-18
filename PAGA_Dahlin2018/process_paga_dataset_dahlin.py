# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import scanpy.api as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
from scipy.stats import mode
from collections import Counter

#you really need anndata 0.6.22, scanpy == 1.4.2

sc.settings.verbosity = 3
sc.set_figure_params(color_map='viridis')
sc.logging.print_versions()
results_file = '/Users/jasoncosgrove/Desktop/dahlin.h5ad'


#this is if you want to send the raw data to file to read into Seurat
#adata = sc.read('/Users/jasoncosgrove/Dropbox/Research/Postdoc_Paris/PAGA_datasets/dahlin18_theis/scanpy_filtered_all_genes_matrix_both.h5ad')
#adata.var_names_make_unique()
#sc.pp.log1p(adata)
#adata.write('/Users/jasoncosgrove/Dropbox/Research/Postdoc_Paris/PAGA_datasets/dahlin18_theis/dahlin_raw.h5ad')



adata = sc.read('/Users/jasoncosgrove/Dropbox (Team Perié)/Jason/Experiments/JCB4_metabolism_bioinformatics_pipeline/datasets/Dahlin/theis_preprocessing/scanpy_processed_matrix_both_no_cc.h5ad')

adata.var_names_make_unique()

#adata.X = None

adata_raw = sc.read('/Users/jasoncosgrove/Dropbox (Team Perié)/Jason/Experiments/JCB4_metabolism_bioinformatics_pipeline/datasets/Dahlin/theis_preprocessing/scanpy_filtered_all_genes_matrix_both.h5ad')
adata_raw.var_names_make_unique()
sc.pp.log1p(adata_raw)
adata.raw = adata_raw
#adataprocessed.write('/Users/jasoncosgrove/Dropbox/Research/Postdoc_Paris/PAGA_datasets/dahlin18_theis/dahlin_processed.h5ad')


sc.pp.neighbors(adata, n_neighbors=7, n_pcs=20)

sc.tl.louvain(adata)

sc.tl.paga(adata)

sc.pl.paga(adata, color=['louvain', 'Procr', 'Gata1', 'Elane'], edge_width_scale=0.2, threshold=0.2)

sc.tl.umap(adata)

adata.obsm['X_umap'][:, 1] = np.clip(adata.obsm['X_umap'][:, 1], a_min=-20, a_max=None)

sc.pl.umap(adata, color=['louvain', 'Procr', 'Gata1', 'Elane'], legend_loc='on data')


sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_diffmap')

sc.tl.louvain(adata, resolution=0.2)

sc.pl.umap(adata, color='louvain', legend_loc='on data')


sc.tl.louvain(adata, resolution=0.06, restrict_to=('louvain', ['4','10','12']))




adata.obs['louvain_R'].cat.categories = [str(i) for i in range(len(adata.obs['louvain_R'].cat.categories))]

adata.obs['clusters'] = adata.obs['louvain_R']

adata.obs['louvain_R'].cat.categories

sc.pl.umap(adata, color=['louvain_R', 'Itga2b', 'Cst3', 'Cd19', 'Dntt'], legend_loc='on data')
sc.pl.umap(adata, color=['louvain_R', 'Procr', 'Gata1', 'Elane', 'Irf8'], legend_loc='on data')
adata.obs['clusters'].cat.categories = [
    '0/MPP', '1/CLP', '2/MPP', '3/E', '4/GMP', '5/Neu', '6/Mono', '7/Baso', '8/CMP', '9/E', '10/MEP', '11/E', '12/Mk',
    '13/Stem', '14/B','15/DC']


sc.tl.paga(adata, groups='clusters')


sc.pl.paga(adata, threshold=0.02, edge_width_scale=1, layout='fr', random_state=7)





#sc.tl.umap(adata, init_pos='paga')

#sc.pl.umap(adata, color='clusters', legend_loc='on data')

#you need to re-run PCA or it will throw an error when you import into R
sc.tl.pca(adata)

#sc.pp.scale(adata)
adata.write(results_file)

#adata.uns['pca']['variance']




adata = sc.read(results_file)
root_cell = int(np.where(adata.raw[:, 'Procr'].X == float(np.max(adata.raw[:, 'Procr'].X)))[0])
print(adata.obs['clusters'][root_cell])
adata.uns['iroot'] = root_cell

sc.tl.dpt(adata)

sc.pl.umap(adata, color=['clusters', 'dpt_pseudotime'], legend_loc='on data')


#gene_names = ['Sell',"Emb","Icam4","Fermt3"]                     # monocyte
gene_names = ["Gata1","Slc14a1","Atpif1","Gstm5","Uqcr11","Abcb10","Gsr","Gpx1","Idh1","Hexa","Hk3","Pkm","Spi1"]                     # monocyte


paths = [('erythrocytes', [13,2,0,8, 10, 3,9, 11]),
         ('monocytes', [13, 2, 0,8, 4,6]),
         ('basophils', [13, 2, 0,8, 4,7]),
         ('neutrophils', [13, 2, 0,8, 4,5]),]


adata.obs['distance'] = adata.obs['dpt_pseudotime'] # this is the distance from the root cell


adata.raw = sc.AnnData(sc.pp.scale(adata.raw.X, copy=True),
                       var=adata.raw.var)



_, axs = pl.subplots(ncols=4, figsize=(12, 6), gridspec_kw={'wspace': 0.05, 'left': 0.12})
pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)
for ipath, (descr, path) in enumerate(paths):
    _, data = sc.pl.paga_path(
        adata, path, gene_names,                         
        show_node_names=False,
        ax=axs[ipath],
        ytick_fontsize=12,
        left_margin=0.15,
        n_avg=400,
        annotations=['distance'],
        show_yticks=True if ipath==0 else False,
        show_colorbar=False,
        color_map='Purples',
        color_maps_annotations={'distance': 'viridis'},
        title='{} path'.format(descr),
        return_data=True,
        show=False)
    data.to_csv('/users/jasoncosgrove/Desktop/paga_path_{}.csv'.format(descr))
pl.savefig('./figures/paga_path_dahlin18.pdf')
pl.show()
