"""
Caleb Ellington, Anush Devadhasan, Daniel Jeong
"""

import numpy as np
import pandas as pd
from causalnex.structure import dynotears


def learn_dbn(timeseries_df, lags):
    W_dy = dynotears.from_pandas_dynamic(timeseries_df, lags)
    W_dy.threshold_till_dag()
    return W_dy


def largest_subgraph(dbn):
    dbn_new = dbn.copy()
    subgraph = dbn_new.get_largest_subgraph()
    all_nodes = dbn_new.nodes()
    remove = set(all_nodes) - set(subgraph)
    dbn_new.remove_nodes_from(remove)
    return dbn_new


def adjacency_matrix(dbn):
    nodes = dbn.nodes
    W = np.zeros((len(nodes), len(nodes)))
    for i, node in enumerate(nodes):
        children = dbn.adj[node]
        for child, value in children.items():
            j = list(nodes).index(child)
            W[i,j] = value['weight']
    return W


def mse(dbn, timeseries_df, lags):
    # make lag sample df
    genes = timeseries_df.columns
    lag_cols = [gene+f"_lag{i}" for i in range(lags+1) for gene in genes]
    lag_rows = []
    for i in range(len(timeseries_df) - lags - 1):
        row = []
        for j in range(lags + 1):
            row += list(timeseries_df.iloc[i + lags])
        lag_rows.append(row)
    lag_df = pd.DataFrame(lag_rows, columns=lag_cols)
    
    W = adjacency_matrix(dbn)
    X = lag_df[list(dbn.nodes)].to_numpy()

    X_prime = X @ W
    mse_val = np.mean((X - X_prime)**2)
    return mse_val