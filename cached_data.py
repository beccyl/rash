####
# Module cached_data: contains all the cached data, and memoized functions to pull data from redis
####
from app import app

import os
import pyarrow.feather
import pandas as pd
from enum import Enum

from flask_caching import Cache

# potentially this could be changed to include:
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
# 'Blood Derived Normal' == 10
# 'Solid Tissue Normal' == 11
# 'Primary Tumor' == 01
# 'Metastatic' == 06
# for now: (Tumor) (Normal)
class sampleType(Enum):
    tumor = 1
    normal = 2

# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
# THCA = Thyroid Carcinoma
# BRCA = Breast Carcinoma
class TCGAStudyType(Enum):
    thyroid = 1
    breast = 2

#### global variables ####
#### TODO: redefine in config file

#tumor_type = "breast"
#app_prefix = "breast_all"
data_dir = "../data/"

CACHE_CONFIG = {
    'CACHE_TYPE': 'redis',
    'CACHE_REDIS_URL': os.environ.get('REDIS_URL', 'redis://localhost:6379'),
    'CACHE_DEFAULT_TIMEOUT': 18000
}
cache = Cache()
cache.init_app(app.server, config=CACHE_CONFIG)
#cache.clear()

@cache.memoize()
def get_sample_list(tumor_type):
    sample_df = pd.read_feather(data_dir + tumor_type + '_tumor_normal_sample.feather')
    sample_dict = sample_df.to_dict('records')
    return sample_dict


@cache.memoize()
def get_N_agg_df_junctions(tumor_type):
    N_agg_df_junctions = pd.read_feather(data_dir + "N_" + tumor_type + '_junctions.feather')
    N_agg_df_junctions.set_index(["sample","gene","chr","start","end"], inplace=True, drop=False)
    N_agg_df_junctions.sort_index(inplace=True)

    return N_agg_df_junctions


@cache.memoize()
def get_N_agg_df_gene_summary(tumor_type):
    N_agg_df_gene_summary = pd.read_feather(data_dir + "N_" + tumor_type + "_gene_summary.feather")
    N_agg_df_gene_summary.set_index(["sample","gene"], inplace=True, drop=False)
    N_agg_df_gene_summary.sort_index(inplace=True)

    return N_agg_df_gene_summary


@cache.memoize()
def get_T_agg_df_junctions(tumor_type):
    T_agg_df_junctions = pd.read_feather(data_dir + "T_" + tumor_type + "_junctions.feather")
    T_agg_df_junctions.set_index(["sample","gene","chr","start","end"], inplace=True, drop=False)
    T_agg_df_junctions.sort_index(inplace=True)

    return T_agg_df_junctions


@cache.memoize()
def get_T_agg_df_gene_summary(tumor_type):
    T_agg_df_gene_summary = pd.read_feather(data_dir + "T_" + tumor_type + "_gene_summary_with_pseudogene.feather")
    T_agg_df_gene_summary.set_index(["sample","gene"], inplace=True, drop=False)
    T_agg_df_gene_summary.sort_index(inplace=True)

    return T_agg_df_gene_summary


@cache.memoize()
def get_T_subset(tumor_type, sample):
    T_agg_df_junctions = get_T_agg_df_junctions(tumor_type)
    T_subset = T_agg_df_junctions.loc[sample,]
    return T_subset


@cache.memoize()
def get_T_summary(tumor_type, sample):
    T_agg_df_gene_summary = get_T_agg_df_gene_summary(tumor_type)
    T_summary = T_agg_df_gene_summary.loc[sample,]
    return T_summary


@cache.memoize()
def get_N_subset(tumor_type, sample):
    N_agg_df_junctions = get_N_agg_df_junctions(tumor_type)
    N_subset = N_agg_df_junctions.loc[sample,]
    return N_subset


@cache.memoize()
def get_N_summary(tumor_type, sample):
    N_agg_df_gene_summary = get_N_agg_df_gene_summary(tumor_type)
    N_summary = N_agg_df_gene_summary.loc[sample,]
    return N_summary


@cache.memoize()
def get_sample_summary_table(tumor_type):
    T_sample_gb = get_T_agg_df_gene_summary(tumor_type).groupby(level=["sample"])
    T_sample_table = pd.concat([T_sample_gb.sum(), (T_sample_gb.size().rename("gene_count")).to_frame()], axis=1).sort_index()
    N_sample_gb = get_N_agg_df_gene_summary(tumor_type).groupby(level=["sample"])
    N_sample_table = pd.concat([N_sample_gb.sum(), (N_sample_gb.size().rename("gene_count")).to_frame()], axis=1).sort_index()

    df_sample_summary_table = pd.merge(T_sample_table, N_sample_table, how="outer", left_index=True, right_index=True, suffixes=("_T","_N")).fillna(0)
    df_sample_summary_table.reset_index(inplace=True, drop=False)

    return df_sample_summary_table


@cache.memoize()
def get_gene_summary_table(tumor_type):
    T_gene_gb = get_T_agg_df_gene_summary(tumor_type).groupby(level=["gene"])
    T_gene_table = pd.concat([T_gene_gb.sum(), (T_gene_gb.size().rename("sample_count")).to_frame()], axis=1).sort_index()
    N_gene_gb = get_N_agg_df_gene_summary(tumor_type).groupby(level=["gene"])
    N_gene_table = pd.concat([N_gene_gb.sum(), (N_gene_gb.size().rename("sample_count")).to_frame()], axis=1).sort_index()

    df_gene_summary_table = pd.merge(T_gene_table, N_gene_table, how="outer", left_index=True, right_index=True, suffixes=("_T","_N")).fillna(0)
    df_gene_summary_table.reset_index(inplace=True, drop=False)

    return df_gene_summary_table

@cache.memoize()
def get_filtered_junctions_count_table(tumor_type, split_read_min_T, split_pair_min_T, junction_count_T, split_read_min_N, split_pair_min_N):
    T_junctions = get_T_agg_df_junctions(tumor_type)[["split_read_count","split_pair_count"]]
    T_junctions = T_junctions[(T_junctions.split_read_count >= split_read_min_T) & (T_junctions.split_pair_count >= split_pair_min_T)]

    T_j_gb = T_junctions.groupby(level=["sample","gene"])

    df_T = (T_j_gb.size().rename("nr_junction_filtered_T")).to_frame()
    df_T = df_T[df_T.nr_junction_filtered_T >= junction_count_T]

    N_junctions = get_N_agg_df_junctions(tumor_type)
    N_junctions = N_junctions[(N_junctions.split_read_count >= split_read_min_N) & (N_junctions.split_pair_count >= split_pair_min_N)]

    N_j_gb = N_junctions.groupby(level=["sample","gene"])

    df_N = (N_j_gb.size().rename("nr_junction_filtered_N")).to_frame()
    df_N = df_N[df_N.nr_junction_filtered_N >= junction_count_T]

    output_table = pd.concat(
        [
            (df_T.groupby(level=["gene"]).size().rename("filtered_sample_count_T")).to_frame(),
            (df_N.groupby(level=["gene"]).size().rename("filtered_sample_count_N")).to_frame()
        ], axis=1).fillna(0)

    return output_table
