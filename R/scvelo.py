import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

adata = sc.read_h5ad('/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/RNAvelocity/my_data.h5ad')

# load loom files for spliced/unspliced matrices for each sample:
ldata1 = scv.read('/home/data/Shared/AQPhan/hiF-T+A549_scRNA-seq_11-16-20/Velocyto_Out/D6LacZ/possorted_genome_bam_7MVG3_D6LacZ.loom', cache=True)
ldata2 = scv.read('/home/data/Shared/AQPhan/hiF-T+A549_scRNA-seq_11-16-20/Velocyto_Out/D6S3/possorted_genome_bam_GHOLO_D6S3.loom', cache=True)
ldata3 = scv.read('/home/data/Shared/AQPhan/hiF-T+A549_scRNA-seq_11-16-20/Velocyto_Out/D9LacZ/possorted_genome_bam_EZ311_D9LacZ.loom', cache=True)
ldata4 = scv.read('/home/data/Shared/AQPhan/hiF-T+A549_scRNA-seq_11-16-20/Velocyto_Out/D9S3/possorted_genome_bam_JLBHK_D9S3.loom', cache=True)
ldata5 = scv.read('/home/data/Shared/AQPhan/hiF-T+A549_scRNA-seq_11-16-20/Velocyto_Out/D12LacZ/possorted_genome_bam_9IE2M_D12LacZ.loom', cache=True)
ldata6 = scv.read('/home/data/Shared/AQPhan/hiF-T+A549_scRNA-seq_11-16-20/Velocyto_Out/D12S3/possorted_genome_bam_3WJ73_D12S3.loom', cache=True)
ldata7 = scv.read('/home/data/Shared/AQPhan/hiF-T+A549_scRNA-seq_11-16-20/Velocyto_Out/D15LacZ/possorted_genome_bam_1QLD5_D15LacZ.loom', cache=True)
ldata8 = scv.read('/home/data/Shared/AQPhan/hiF-T+A549_scRNA-seq_11-16-20/Velocyto_Out/D15S3/possorted_genome_bam_8JWL0_D15S3.loom', cache=True)

# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_61' for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_62' for bc in barcodes]
ldata2.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_91' for bc in barcodes]
ldata3.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata4.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_92' for bc in barcodes]
ldata4.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata5.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_121' for bc in barcodes]
ldata5.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_122' for bc in barcodes]
ldata3.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_151' for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_152' for bc in barcodes]
ldata2.obs.index = barcodes

# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()
# concatenate the three loom
ldata = ldata1.concatenate([ldata2, ldata3])
# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)
# plot umap to check
sc.pl.umap(adata, color='celltype', frameon=False, legend_loc='on data', title='', save='_celltypes.pdf')
