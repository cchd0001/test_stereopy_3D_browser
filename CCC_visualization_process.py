# -*- coding: utf-8 -*-
# @Time    : 2023/6/6 13:41
# @Author  : liuxiaobin
# @File    : CCC_visualization_process.py
# @Version：V 0.1
# @desc :

import pandas as pd
import numpy as np
import anndata
import sys

def get_cell_type_1(adata):
    """
    Return list of all possible celltype1.
    """
    celltype1 = list(set(adata.uns['ccc_data']['celltype1']))
    celltype1.sort()
    return celltype1


def get_cell_type_2(adata, celltype1):
    """
    Given a celltype1, return a list of all possible celltype2.
    """
    data = adata.uns['ccc_data']
    data1 = data[data['celltype1'] == celltype1]
    celltype2 = list(set(data1['celltype2']))
    celltype2.sort()
    return celltype2


def get_ligand(adata, celltype1, celltype2):
    """
    Given a celltype1 and a celltype2, return a list of all possible ligands.
    """
    data = adata.uns['ccc_data']
    data12 = data[(data['celltype1'] == celltype1) & (data['celltype2'] == celltype2)]
    ligand = list(set(data12['ligand']))
    ligand.sort()
    return ligand


def get_receptor(adata, celltype1, celltype2, ligand):
    """
    Given a celltype1, a celltype2 and a ligand, return a list of all possible receptors.
    """
    data = adata.uns['ccc_data']
    data12 = data[(data['celltype1'] == celltype1) & (data['celltype2'] == celltype2) & (data['ligand'] == ligand)]
    receptor = list(set(data12['receptor']))
    receptor.sort()
    return receptor


def get_data_for_CCC_visulization(adata, celltype1, celltype2, ligand, receptor, celltype_col='celltype'):
    """
    Given celltype pairs and lr pairs, return expression data for visualization.
    """
    adata_ligand = adata[adata.obs[celltype_col] == celltype1, adata.var.index.str.lower() == ligand.lower()]
    adata_receptor = adata[adata.obs[celltype_col] == celltype2, adata.var.index.str.lower() == receptor.lower()]
    return adata_ligand, adata_receptor
