"""
Caleb Ellington, Anush Devadhasan, Daniel Jeong
"""

import csv
import gzip
import os
import os.path as osp

import numpy as np
import scipy.io
import scanpy as sc


def load_raw_spots(mat_dir='./data/filtered_feature_bc_matrix'):
    """
    todo: docstring
    """
    mat_path = osp.join(mat_dir, 'matrix.mtx.gz')
    features_path = osp.join(mat_dir, 'features.tsv.gz')
    barcodes_path = osp.join(mat_dir, 'barcodes.tsv.gz')
    
    mat = scipy.io.mmread(mat_path).todense() # (36601, 4325)
    feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, 'rt'), delimiter='\t')] # Ensembl IDs
    gene_names = [row[1] for row in csv.reader(gzip.open(features_path, 'rt'), delimiter='\t')] # Gene Names
    barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, 'rt'), delimiter='\t')] # Spot Barcodes
    return mat, gene_names, barcodes


def load_processed_spots(sample_id='Parent_Visium_Human_BreastCancer', n_genes=30):
    """
    todo: docstring
    """
    # Load the data
    adata = sc.datasets.visium_sge(sample_id=sample_id)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-") # Prefix MT- marks mitochondrial genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True) # Calculates some Quality Control metrics
    
    # Filter and Normalize the data
    sc.pp.filter_cells(adata, min_counts=5000)
    sc.pp.filter_cells(adata, max_counts=35000)
    adata = adata[adata.obs["pct_counts_mt"] < 20] # Only take spots without overly expressed mitochondrial genes
    sc.pp.filter_genes(adata, min_cells=10) # Filter out genes with less than total counts of 10
    sc.pp.normalize_total(adata, inplace=True) # Normalizes total count to median of total counts for observations (cells) before normalization
    sc.pp.log1p(adata) # log-transform the expression data i.e. for X, calculate log(X+1)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=n_genes)  # Take <n_genes> most interesting genes
    
    # Format the data and return
    mat = adata.X[:, adata.var['highly_variable']].todense()
    gene_names = list(adata.var['gene_ids'][adata.var['highly_variable']].index)
    barcodes = list(adata.obs.index)
    return mat, gene_names, barcodes
    
    
def load_positions(barcodes, positions_path='./data/spatial/tissue_positions_list.csv'):
    """
    todo: docstring
    """
    positions_dict = {line[0]:np.array(line)[1:].astype(np.int32) for line in csv.reader(open(positions_path, 'r'))}
    positions = np.array([positions_dict[barcode] for barcode in barcodes]) # Spot positions (row, col, px_row, px_col)
    return positions
 