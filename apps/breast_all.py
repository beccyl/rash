#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
import os
import re
import sys
import pyarrow.feather

import pandas as pd
#import numpy as np
from scipy import stats
from enum import Enum

import plotly.graph_objs as go
from plotly.subplots import make_subplots

import plotly.express as px
import plotly.figure_factory as ff

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table
from dash.dependencies import Input, Output
from app import app
from table_bars import data_bars

from flask_caching import Cache

#import cached_data
from cached_data import *

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class filterLevel(Enum):
    gene = 1
    junction = 2

#### global variables ####
#### TODO: redefine in config file
tumor_type = "breast"
app_prefix = tumor_type + "_all_"

# default values
filter_level_default_T = filterLevel.junction.name
filter_level_default_N = filterLevel.junction.name
split_read_min_default = 1
split_pair_min_default = 0
junctions_default = 2
coverage_default = None  ## TODO: change this

#sample_dict = get_sample_list(tumor_type)

############## load initial data ##### TODO: move this!!!
def safe_eval(x):
    if (isinstance(x, int)):
        return float(x)
    try:
        return eval(x)
    except ZeroDivisionError:
        return float("inf")

#######################################
## Classes                           ##
#######################################

# read json file
# create dataframe from json file

class Sample:
    def __init__(self, sample_name, logger):
        self.logger = logger
        self.sample_name = sample_name
            #self.read_json() # can raise IOError or FileNotFoundException
        self.basename = os.path.splitext(os.path.splitext(os.path.basename(self.sample_name))[0])[0]
        logger.info("basename:[%s]", self.basename)
        logger.info("filename:[%s]", self.sample_name)

    def getJunctions(self):
        return self.df_junctions

    def setBasename(self, basename):
        self.basename = basename

    def setFilename(self, sample_name):
        self.sample_name = sample_name

    def setFrames(self, df, df_junctions):
        self.df_junctions = df_junctions
        self.df = df

    def getGeneNames(self):
        return self.df["gene"]

    def getGeneJunctions(self,gname=None):
        if gname is not None:
            return self.df[self.df["gene"] == gname][["gene","junction_count"]]
        return self.df[["gene","junction_count"]]

    def reportJunctionComparison(
            self,
            other=None,
            level_T=filterLevel.gene,
            #coverage_N=None,  #TODO: filter normal by this (?)
            level_N=filterLevel.gene,
            split_pair_min_N=0,
            split_read_min_N=0,
            #nr_junct_N=1,     #TODO: filter normal by this (?)
            coverage_T=None,
            split_pair_min_T=0,
            split_read_min_T=0,
            nr_junct_T=1,
            pseudo_df=None,
            value_counts=None):

        dfSelf = self.df  # default
        df_junctions = self.df_junctions #default
        dfForCSV = dfSelf

        self_index_names = dfSelf.index.names
        self_junctions_index_names = self.df_junctions.index.names

        ### TODO: Simplify
        if other is not None and other.df is not None and other.df_junctions is not None:
            other_index_names = other.df.index.names
            other_junctions_index_names = other.df_junctions.index.names

            if set(dfSelf.index.names) != set(other.df.index.names):
                try:
                    dfSelf = dfSelf.reset_index(level=list(set(dfSelf.index.names)-set(other.df.index.names)),drop=True)
                    logger.info("dfSelf.index.names %s", dfSelf.index.names)
                    logger.info("other.df.index.names %s", other.df.index.names)
                except:
                    logger.info("BUG: TODO: FIXME")

            if set(df_junctions.index.names) != set(other.df_junctions.index.names):
                try:
                    df_junctions = df_junctions.reset_index(level=list(set(df_junctions.index.names)-set(other.df_junctions.index.names)), drop=True)
                    logger.info("df_junctions.index.names %s", df_junctions.index.names)
                    logger.info("other.df_junctions.index.names %s", other.df_junctions.index.names)
                except:
                    logger.info("BUG: TODO: FIXME")


        if level_N == filterLevel.gene:
            if other is not None and other.df is not None:
                index = other.df[(other.df.split_read_count >= split_read_min_N) & (other.df.split_pair_count >= split_pair_min_N)].index

                dfSelf = dfSelf[~dfSelf.index.isin(index)]

        if level_T == filterLevel.gene:
            dfSelf = dfSelf[(dfSelf.split_read_count >= split_read_min_T) & (dfSelf.split_pair_count >= split_pair_min_T) & (dfSelf.junction_count >= nr_junct_T)]


            #TODO: Simplify
            self_junctions_index_names = df_junctions.index.names
            if set(dfSelf.index.names) != set(df_junctions.index.names):
                logger.info("dfSelf.index.names %s", dfSelf.index.names)
                logger.info("df_junctions.index.names %s", df_junctions.index.names)
                try:
                    df_junctions = df_junctions.reset_index(level=list(set(dfSelf.index.names)-set(df_junctions.index.names)), drop=True)
                except:
                    logger.info("BUG: TODO: FIXME")
                df_junctions.set_index(self_junctions_index_names, drop=False)
            #df_junctions = df_junctions.reset_index(level=["chr","start","end"], drop=True)
            #df_junctions = df_junctions.loc[dfSelf.index,:]
            #df_junctions = df_junctions.set_index(["chr","start","end"], append=True, drop=False)

        if level_N == filterLevel.junction:
            if other is not None and other.df_junctions is not None:
                index = other.df_junctions[(other.df_junctions.split_read_count >= split_read_min_N) & (other.df_junctions.split_pair_count >= split_pair_min_N)].index

                df_junctions = df_junctions[~df_junctions.index.isin(index)]

        if level_T == filterLevel.junction:
            df_junctions = df_junctions[(df_junctions.split_read_count >= split_read_min_T) & (df_junctions.split_pair_count >= split_pair_min_T)]

        if coverage_T is not None:
            df_junctions = df_junctions[df_junctions.ds_coverage.apply(lambda x: safe_eval(x) < coverage_T) & df_junctions.us_coverage.apply(lambda x: safe_eval(x) < coverage_T)]

        ### if the index was dropped -- reinstate it ...
        # TODO: Simplify
        if set(self_index_names) != set(dfSelf.index.names):
            dfSelf = dfSelf.set_index(list(set(self_index_names)-set(dfSelf.index.names)), append=True, drop=False)
            dfSelf = dfSelf.reorder_levels(self_index_names, axis="index")
            dfSelf.sort_index(inplace=True)

        if set(self_junctions_index_names) != set(df_junctions.index.names):
            df_junctions = df_junctions.set_index(list(set(self_junctions_index_names)-set(df_junctions.index.names)), append=True, drop=False)
            df_junctions = df_junctions.reorder_levels(self_junctions_index_names, axis="index")
            df_junctions.sort_index(inplace=True)

        ## group junctions by levels in gene summary (ie. ["sample", "gene"], or ["gene"])
        junction_counts = df_junctions.groupby(level=dfSelf.index.names).size()
        junction_counts.rename("nr_junct", inplace=True)

        ## ensure that the number of junctions meets criteria
        dfForCSV = dfSelf[dfSelf.index.isin(junction_counts[junction_counts >= nr_junct_T].index)]
        dfForCSV = pd.merge(dfForCSV, junction_counts, how="left", left_index=True, right_index=True).fillna(0)

        if logger.isEnabledFor(logging.INFO):
            logger.info("Filter all T=%s N=%s: [split_pair_min_T=%s], [split_read_min_T=%s], [nr_junct_T=%s], [coverage_T=%s], [split_pair_min_N=%s], [split_read_min_N=%s]"#, [nr_junct_N=%s], [coverage_N=%s]"
                , level_T, level_N, str(split_pair_min_T), str(split_read_min_T), str(nr_junct_T), str(coverage_T)
                , str(split_pair_min_N), str(split_read_min_N))#, str(nr_junct_N), str(coverage_N))

        #####
        ## create a copy rather than a view  -- this is important!! if not a copy, will bug
        #####
        dfForCSV = dfForCSV.assign(id=dfForCSV.index.to_flat_index())
        if pseudo_df is not None:
            dfForCSV.insert(0, "pseudogenes","") ## empty -- to override
        if value_counts is not None:
            dfForCSV.insert(0, "T_value_count",0) ## empty -- to override
            dfForCSV.insert(0, "N_value_count",0) ## empty -- to override

        if pseudo_df is not None:
            grouped = pseudo_df.groupby(by="parent_gene_name")

            for index, row in dfForCSV.iterrows():
                if value_counts is not None:
                    dfForCSV.loc[index, "T_value_count"] = value_counts.loc[row.gene, "T_counts"]
                    dfForCSV.loc[index, "N_value_count"] = value_counts.loc[row.gene, "N_counts"]
                try:
                    values = grouped.get_group(row.gene).pseudogene_gene_name.to_string(header=False,index=False)
                    values = ",".join(values.split())
                    dfForCSV.loc[index, "pseudogenes"] = values
                except KeyError:
                    pass

        dfForCSV = pd.merge(dfForCSV.reset_index(drop=True), get_gene_summary_table(tumor_type)[["gene","sample_count_T","sample_count_N"]], how="left", left_on="gene", right_on="gene").fillna(0)

        dfForCSV = pd.merge(dfForCSV.reset_index(drop=True), get_N_agg_df_gene_summary(tumor_type).reset_index(drop=True), how="left", left_on=["sample","gene"], right_on=["sample","gene"], suffixes=("_T","_N")).fillna(0)

        return dfForCSV

########################
#### END CLASS DEFN ####
########################

@app.callback(
    Output(app_prefix+'filter_level_T', 'options'),
    [
        Input(app_prefix+'filter_status_selector', 'value')
    ]
)
def set_filter_level_T_options(filter_status_selector):
    if filter_status_selector == "all":
        return [{"label": filterLevel.junction.name, "value": filterLevel.junction.name}]
    else:
        return [
            {"label": filterLevel.gene.name, "value": filterLevel.gene.name},
            {"label": filterLevel.junction.name, "value": filterLevel.junction.name}
        ]


@app.callback(
    Output(app_prefix+'filter_level_N', 'options'),
    [
        Input(app_prefix+'filter_status_selector', 'value')
    ]
)
def set_filter_level_N_options(filter_status_selector):
    if filter_status_selector == "all":
        return [{"label": filterLevel.junction.name, "value": filterLevel.junction.name}]
    else:
        return [
            {"label": filterLevel.gene.name, "value": filterLevel.gene.name},
            {"label": filterLevel.junction.name, "value": filterLevel.junction.name}
        ]


@app.callback(
    Output(app_prefix + 'filter_level_T', 'value'),
    [
        Input(app_prefix + 'filter_level_T', 'options')
    ]
)
def set_filter_level_T_value(available_options):
    return available_options[len(available_options)-1]["value"]


@app.callback(
    Output(app_prefix + 'filter_level_N', 'value'),
    [
        Input(app_prefix + 'filter_level_N', 'options')
    ]
)
def set_filter_level_N_value(available_options):
    return available_options[len(available_options)-1]["value"]
    #return available_options[0]["value"]


@app.callback(
    Output(app_prefix + 'gene_summary_table', 'page_current'),
    [
        #Input(app_prefix + 'sample_dropdown', 'value'),
        Input(app_prefix + 'filter_status_selector', 'value'),
        Input(app_prefix + 'filter_level_T', 'value'),
        Input(app_prefix + 'split_read_min_T', 'value'),
        Input(app_prefix + 'split_pair_min_T', 'value'),
        Input(app_prefix + 'junctions_T', 'value'),
        Input(app_prefix + 'coverage_T', 'value'),
        Input(app_prefix + 'coverage_check_T', 'value'),
        Input(app_prefix + 'filter_level_N', 'value'),
        Input(app_prefix + 'split_read_min_N', 'value'),
        Input(app_prefix + 'split_pair_min_N', 'value'),
    ]
)
def update_page_current(
    #sample_dropdown,
    filter_status_selector,
    filter_level_T,
    split_read_min_T,
    split_pair_min_T,
    junctions_T,
    coverage_T,
    coverage_check_T,
    filter_level_N,
    split_read_min_N,
    split_pair_min_N):#,
    return 0


@app.callback(
    Output(app_prefix + 'gene_summary_table', 'derived_virtual_selected_row_ids'),
    [
        #Input(app_prefix + 'sample_dropdown', 'value'),
        Input(app_prefix + 'filter_status_selector', 'value'),
        Input(app_prefix + 'filter_level_T', 'value'),
        Input(app_prefix + 'split_read_min_T', 'value'),
        Input(app_prefix + 'split_pair_min_T', 'value'),
        Input(app_prefix + 'junctions_T', 'value'),
        Input(app_prefix + 'coverage_T', 'value'),
        Input(app_prefix + 'coverage_check_T', 'value'),
        Input(app_prefix + 'filter_level_N', 'value'),
        Input(app_prefix + 'split_read_min_N', 'value'),
        Input(app_prefix + 'split_pair_min_N', 'value'),
    ]
)
def update_selected_virtual_row_ids(
    #sample_dropdown,
    filter_status_selector,
    filter_level_T,
    split_read_min_T,
    split_pair_min_T,
    junctions_T,
    coverage_T,
    coverage_check_T,
    filter_level_N,
    split_read_min_N,
    split_pair_min_N):#,
    return []


@app.callback(
    Output(app_prefix + 'gene_summary_table', 'selected_rows'),
    [
        #Input(app_prefix + 'gene_summary_table', 'selected_row_ids'),
        #Input(app_prefix + 'sample_dropdown', 'value'),
        Input(app_prefix + 'filter_status_selector', 'value'),
        Input(app_prefix + 'filter_level_T', 'value'),
        Input(app_prefix + 'split_read_min_T', 'value'),
        Input(app_prefix + 'split_pair_min_T', 'value'),
        Input(app_prefix + 'junctions_T', 'value'),
        Input(app_prefix + 'coverage_T', 'value'),
        Input(app_prefix + 'coverage_check_T', 'value'),
        Input(app_prefix + 'filter_level_N', 'value'),
        Input(app_prefix + 'split_read_min_N', 'value'),
        Input(app_prefix + 'split_pair_min_N', 'value'),
        #Input(app_prefix + 'junctions_N', 'value'),
        #Input(app_prefix + 'coverage_N', 'value'),
        #Input(app_prefix + 'coverage_check_N', 'value'),
    ]
)
def update_selected_row_ids(
    #sample_dropdown,
    filter_status_selector,
    filter_level_T,
    split_read_min_T,
    split_pair_min_T,
    junctions_T,
    coverage_T,
    coverage_check_T,
    filter_level_N,
    split_read_min_N,
    split_pair_min_N):#,
    #junctions_N,
    #coverage_N,
    #coverage_check_N):
    return []


@app.callback(
    Output(app_prefix + 'gene_summary_table', 'data'),
    [
        #Input(app_prefix + 'sample_selector', 'value'),
        #Input(app_prefix + 'sample_dropdown', 'value'),
        Input(app_prefix + 'filter_status_selector', 'value'),
        Input(app_prefix + 'filter_level_T', 'value'),
        Input(app_prefix + 'split_read_min_T', 'value'),
        Input(app_prefix + 'split_pair_min_T', 'value'),
        Input(app_prefix + 'junctions_T', 'value'),
        Input(app_prefix + 'coverage_T', 'value'),
        Input(app_prefix + 'coverage_check_T', 'value'),
        Input(app_prefix + 'filter_level_N', 'value'),
        Input(app_prefix + 'split_read_min_N', 'value'),
        Input(app_prefix + 'split_pair_min_N', 'value'),
        #Input('junctions_N', 'value'),
        #Input('coverage_N', 'value'),
        #Input('coverage_check_N', 'value'),
    ]
)
@cache.memoize() # in seconds #TODO:
def update_table_data(
        #sample_selector,
        #sample,
        filter_status_selector,
        filter_level_T,
        split_read_min_T,
        split_pair_min_T,
        junctions_T,
        coverage_T,
        coverage_check_T,
        filter_level_N,
        split_read_min_N,
        split_pair_min_N):#,
        #junctions_N,
        #coverage_N,
        #coverage_check_N):

    tumor = None
    normal = None
    #T_subset = get_T_subset(tumor_type, sample)
    #dfT = get_T_summary(tumor_type, sample)
    #N_subset = get_N_subset(tumor_type, sample)
    #dfN = get_N_summary(tumor_type, sample)

    tumor = Sample("all", logger)  ## set the basename and filename to sample
    tumor.setFrames(get_T_agg_df_gene_summary(tumor_type), get_T_agg_df_junctions(tumor_type))

    normal = Sample("all", logger)
    normal.setFrames(get_N_agg_df_gene_summary(tumor_type), get_N_agg_df_junctions(tumor_type))


    if filter_status_selector == "None":
        normal = None
    elif filter_status_selector == "sample":
        normal = normal
    elif filter_status_selector == "all":
        normal_aggregation = Sample("all", logger)
        normal_aggregation.setFrames(
            get_N_agg_df_gene_summary(tumor_type).reset_index(level="sample", drop=True),
            get_N_agg_df_junctions(tumor_type).reset_index(level="sample", drop=True),
        )
        normal = normal_aggregation

    logger.info("filter_status_selector: [%s]", filter_status_selector)

    # TODO: use safe_eval function (?)
    if not coverage_check_T:
        coverage_T = None
    else:
        try:
            coverage_T = float(coverage_T)
        except TypeError:
            logger.warn("coverage_T: [%s] is not a float", coverage_T)
            coverage_T = None


    df_for_table = tumor.reportJunctionComparison(
        normal,
        level_T=filterLevel[filter_level_T],
        split_read_min_T=int(split_read_min_T),
        split_pair_min_T=int(split_pair_min_T),
        level_N=filterLevel[filter_level_N],
        split_read_min_N=int(split_read_min_N),
        split_pair_min_N=int(split_pair_min_N),
        coverage_T=coverage_T,
        nr_junct_T=int(junctions_T)
    )

    df_agg_count_filtered = get_filtered_junctions_count_table(tumor_type, split_read_min_T, split_pair_min_T, junctions_T, split_read_min_N, split_pair_min_N)

    df_output = pd.merge(df_for_table, df_agg_count_filtered, how="left", left_on="gene", right_index=True).fillna(0)

    #return df_for_table.to_dict('records')
    return df_output.to_dict('records')


@app.callback(
    Output(app_prefix + 'junction_table','page_current'),
    [
        Input(app_prefix + 'gene_summary_table','derived_virtual_selected_row_ids'),
    ]
)
def update_junction_table_page_current(selected_row_ids):
    return 0


@app.callback(
    Output(app_prefix + 'junction_table', 'style_data_conditional'),
    [
        Input(app_prefix + 'split_read_min_T','value'),
        Input(app_prefix + 'split_pair_min_T','value')
    ]
)
def update_style_junction_table(
    split_read_min,
    split_pair_min,
):
    return ([{
        'if': {
            'filter_query': '{{split_read_count}} >= {} && {{split_pair_count}} >= {}'.format(split_read_min,split_pair_min)
        },
        'backgroundColor': 'dodgerblue',
        'color': 'white',
        'fontWeight':'bold',
    }] + [{
        'if': {
            'column_id': 'split_read_count',
            'filter_query': '{{split_read_count}} >= {}'.format(split_read_min)
        },
    #    'backgroundColor': 'dodgerblue',
    #    'color': 'white',
        'fontWeight':'bold',
    }] + [{
        'if': {
            'column_id': 'split_pair_count',
            'filter_query': '{{split_pair_count}} >= {}'.format(split_pair_min)
        },
    #    'backgroundColor': 'dodgerblue',
    #    'color': 'white',
        'fontWeight':'bold',
    }])


@app.callback(
    Output(app_prefix + 'junction_table','data'),
    [
        #Input(app_prefix + 'sample_dropdown','value'),
        Input(app_prefix + 'gene_summary_table','derived_virtual_selected_row_ids'),
        #Input(app_prefix + 'filter_status_selector', 'value'),
        #Input(app_prefix + 'split_read_min', 'value'),
        #Input(app_prefix + 'split_pair_min', 'value'),
        #Input(app_prefix + 'junctions', 'value'),
        #Input(app_prefix + 'coverage', 'value')
    ]
)
@cache.memoize() # in seconds #TODO:
def update_junction_table(
    #sample,
    selected_row_ids): #, filter_status_selector, split_read_min, split_pair_min,     junctions, coverage):
    #selected_id_set = set(selected_row_ids or [])

    if not selected_row_ids:
        return []

    selected_id_set = list(map(tuple, selected_row_ids))
    #print(selected_id_set)
    #print(get_T_agg_df_junctions(tumor_type).index.names)

    T_subset = (get_T_agg_df_junctions(tumor_type).reset_index(level=["chr","start","end"],drop=True)).loc[selected_id_set , ]
    #get_T_subset(tumor_type, sample)
    #T_subset_selected = T_subset[T_subset.index.isin(selected_id_set,level="gene")]

    return T_subset.to_dict('records')


@app.callback(
    Output(app_prefix + "coverage_N", "disabled"),
    [Input(app_prefix + 'coverage_check_N', 'value')]
)
def disableCoverage(coverage_check_N):
    return coverage_check_N==""  # if coverage_check is "" (None), then disable the coverage


@app.callback(
    Output(app_prefix + "coverage_T", "disabled"),
    [Input(app_prefix + 'coverage_check_T', 'value')]
)
def disableCoverage(coverage_check_T):
    return coverage_check_T==""  # if coverage_check is "" (None), then disable the coverage


layout = html.Div(
    [
        html.Div(
            [
                html.Div(
                    [
                        html.Img(
                            src=app.get_asset_url("logo.svg"),
                            id="pmac-image",
                            style={
                                "height": "60px",
                                "width": "auto",
                                "margin-bottom": "25px",
                            },
                        )
                    ],
                    className="one-third column",
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.H3(
                                    "Retrogene Analysis",
                                    style={"margin-bottom": "0px"},
                                ),
                                html.H5(
                                    tumor_type.title() + " summary", style={"margin-top": "0px"}
                                ),
                            ]
                        )
                    ],
                    className="one-half column",
                    id="title",
                ),
                html.Div(
                    [
                        html.A(
                            html.Button("View on GitHub", id="learn-more-button"),
                            href="https://github.com/beccyl/",
                        ),
                        #html.Img(src=app.get_asset_url("GitHub-Mark-Light-64px.png")),
                    ],
                    className="one-third column",
                    id="button",
                ),
            ],
            id="header",
            className="row flex-display",
            style={"margin-bottom": "25px"},
        ),
        html.Div(
            [
                html.Div (
                    [
                        html.P("Filter tumor data against:", className="control_label",),
                        dcc.RadioItems(
                            id=app_prefix + "filter_status_selector",
                            options=[
                                {"label": "None", "value": "None"},
                                {"label": "matched normal", "value": "normal"},
                                {"label": "all normals", "value": "all"}
                            ],
                            value="normal",
                            className="dcc_control",
                            #labelStyle={"display":"inline-block"}
                        ),
                        dcc.Tabs(
                            id=app_prefix + "tabs-filter",
                            value="tab-tumor",
                            parent_className="custom-tabs",
                            className="custom-tabs-container",
                            children=[
                            dcc.Tab(
                                label="Tumor filters",
                                value="tab-tumor",
                                children=[
                                    html.P("Filter level:", className="control_label"),
                                    dcc.RadioItems(
                                        id=app_prefix + "filter_level_T",
                                        options=[
                                            {"label": filterLevel.gene.name, "value": filterLevel.gene.name},
                                            {"label": filterLevel.junction.name, "value": filterLevel.junction.name}
                                        ],
                                        value=filter_level_default_T,
                                        className="dcc_control",
                                        labelStyle={"display":"inline-block"}
                                    ),
                                    html.P("split read minimum:", className="control_label"),
                                    dcc.Slider(
                                        id=app_prefix + "split_read_min_T",
                                        min=0,
                                        max=10,
                                        step=1,
                                        dots=True,
                                        marks={
                                            0: "0",
                                            1: "1",
                                            2: "2",
                                            5: "5",
                                            10: "10",
                                        },
                                        value=split_read_min_default,
                                        className="dcc_control",
                                    ), # split_read_min - default 1, 0..N
                                    html.P("split pair minimum:", className="control_label"),
                                    dcc.Slider(
                                        id=app_prefix + "split_pair_min_T",
                                        min=0,
                                        max=10,
                                        step=1,
                                        dots=True,
                                        marks={
                                            0: "0",
                                            1: "1",
                                            2: "2",
                                            5: "5",
                                            10: "10",
                                        },
                                        value=split_pair_min_default,
                                        className="dcc_control",
                                    ), # split_pair_min - default 0, 0..N
                                    html.P("junctions:", className="control_label"),
                                    dcc.Slider(
                                        id=app_prefix + "junctions_T",
                                        min=0,
                                        max=10,
                                        step=None,
                                        marks={
                                            1: "1",
                                            2: "2",
                                            3: "",
                                            4: "",
                                            5: "5",
                                            6: "",
                                            7: "",
                                            8: "",
                                            9: "",
                                            10: "10",
                                        },
                                        value=junctions_default,
                                        className="dcc_control",
                                    ), # junctions - default 1, 1..N
                                    html.P("Coverage drop at junction:", className="control_label"),
                                    dcc.RadioItems(
                                        id=app_prefix + "coverage_check_T",
                                        options=[
                                            {"label": "None", "value": ""},
                                            {"label": "Value", "value": "value"},
                                        ],
                                        value="",
                                        className="dcc_control",
                                        labelStyle={"display":"inline-block"}
                                    ),
                                    dcc.Slider(
                                        id=app_prefix + "coverage_T",
                                        min=0.5,
                                        max=1,
                                        step=None,
                                        marks={
                                            0.5: "0.5",
                                            0.6: "0.6",
                                            0.7: "0.7",
                                            0.8: "0.8",
                                            0.9: "0.9",
                                            1: "1",
                                        },
                                        disabled=True,
                                        value=0.8,
                                        className="dcc_control",
                                    )
                                ]
                            ),
                            dcc.Tab(
                                label="Normal filters",
                                value="tab-normal",
                                children=[
                                    html.P("Filter out normals with criteria:", className="control_label",),
                                    dcc.RadioItems(
                                        id=app_prefix + "filter_level_N",
                                        options=[
                                            {"label": filterLevel.gene.name, "value": filterLevel.gene.name},
                                            {"label": filterLevel.junction.name, "value": filterLevel.junction.name}
                                        ],
                                        value=filter_level_default_N,
                                        className="dcc_control",
                                        labelStyle={"display":"inline-block"}
                                    ),
                                    html.P("split read minimum:", className="control_label"),
                                    dcc.Slider(
                                        id=app_prefix + "split_read_min_N",
                                        min=0,
                                        max=10,
                                        step=1,
                                        dots=True,
                                        marks={
                                            0: "0",
                                            1: "1",
                                            2: "2",
                                            5: "5",
                                            10: "10",
                                        },
                                        value=split_read_min_default,
                                        className="dcc_control",
                                    ), # split_read_min - default 1, 0..N
                                    html.P("split pair minimum:", className="control_label"),
                                    dcc.Slider(
                                        id=app_prefix + "split_pair_min_N",
                                        min=0,
                                        max=10,
                                        step=1,
                                        dots=True,
                                        marks={
                                            0: "0",
                                            1: "1",
                                            2: "2",
                                            5: "5",
                                            10: "10",
                                        },
                                        value=split_pair_min_default,
                                        className="dcc_control",
                                    ), # split_pair_min - default 0, 0..N
                                ]
                            ),
                        ]),
                    ],
                    className="pretty_container three columns",
                    id=app_prefix + "cross-filter-options",
                ),
                html.Div(
                    [
                        html.H6("Genes"),
                        dash_table.DataTable(
                            id=app_prefix + 'gene_summary_table',
                            columns=[
                                dict(id="sample", name=["", "sample"], type="text"),#, presentation="markdown"),
                                dict(id="gene", name=["", "gene"], type="text"),
                                dict(id="split_read_count_T", name=["tumor counts (total)", "split_read_count_T"], type="numeric"),
                                dict(id="split_pair_count_T", name=["tumor counts (total)", "split_pair_count_T"], type="numeric"),
                                #dict(id="junction_count_T", name=["tumor counts (total)", "junction_count_T"], type="numeric"),
                                dict(id="nr_junct", name=["(filtered)", "junctions_T"], type="numeric"),
                                dict(id="split_read_count_N", name=["normal counts (total)", "split_read_count_N"], type="numeric"),
                                dict(id="split_pair_count_N", name=["normal counts (total)", "split_pair_count_N"], type="numeric"),
                                #dict(id="junction_count_N", name=["normal counts (total)", "junction_count_N"], type="numeric"),
                                dict(id="sample_count_T", name=["aggregation", "sample_count_T"], type="numeric"),
                                dict(id="sample_count_N", name=["aggregation", "sample_count_N"], type="numeric"),
                                dict(id="filtered_sample_count_T", name=["filtered", "sample_count_T"], type="numeric"),
                                dict(id="filtered_sample_count_N", name=["filtered", "sample_count_N"], type="numeric"),
                                dict(id="pseudogenes", name=["", "pseudogenes"], type="text"),
                            ],
                            page_size=10,
                            fixed_rows={'headers':True},
                            merge_duplicate_headers=True,
                            editable=False,
                            sort_action="native",
                            #sort_mode="multi",
                            filter_action="native",
                            #row_deletable=False,
                            row_selectable="multi",
                            selected_rows=[],
                            cell_selectable=True,
                            #export_columns="all",
                            export_format="csv",
                            export_headers="names",
                            style_data={
                                'width': '100px', 'minWidth': '100px', 'maxWidth': '100px',
                                'overflow': 'hidden',
                                'textOverflow': 'ellipsis',
                            },
                            style_data_conditional=(
                                [
                                    {
                                        'if': {
                                            'filter_query': '{pseudogenes} is blank',
                                            'column_id': 'pseudogenes',
                                        },
                                        'backgroundColor': 'dodgerblue',
                                        'color': 'white'
                                    }
                                ] +
                                [
                                    {
                                        'if': {
                                            'filter_query': '{split_read_count_T} = 0',
                                            'column_id': 'split_read_count_T',
                                        },
                                        'color': 'tomato',
                                        'fontWeight': 'bold'
                                    }
                                ] +
                                [
                                    {
                                        'if': {
                                            'filter_query': '{split_pair_count_T} = 0',
                                            'column_id': 'split_pair_count_T',
                                        },
                                        'color': 'tomato',
                                        'fontWeight': 'bold'
                                    }
                                ] +
                                [
                                    {
                                        'if': {
                                            'filter_query': '{junction_count_N} > 0',
                                            'column_id': 'junction_count_N',
                                        },
                                        'color': 'tomato',
                                        'fontWeight': 'bold'
                                    }
                                ] +
                                [
                                    {
                                        'if': {
                                            'filter_query': '{split_read_count_N} > 0',
                                            'column_id': 'split_read_count_N',
                                        },
                                        'color': 'tomato',
                                        'fontWeight': 'bold'
                                    }
                                ] +
                                [
                                    {
                                        'if': {
                                            'filter_query': '{split_pair_count_N} > 0',
                                            'column_id': 'split_pair_count_N',
                                        },
                                        'color': 'tomato',
                                        'fontWeight': 'bold'
                                    }
                                ] +
                                data_bars(get_gene_summary_table(tumor_type),'sample_count_T') +
                                data_bars(get_gene_summary_table(tumor_type),'sample_count_N') +
                                data_bars(get_filtered_junctions_count_table(tumor_type, split_read_min_default, split_pair_min_default, junctions_default, split_read_min_default, split_pair_min_default) ,'filtered_sample_count_T') +
                                data_bars(get_filtered_junctions_count_table(tumor_type, split_read_min_default, split_pair_min_default, junctions_default, split_read_min_default, split_pair_min_default) ,'filtered_sample_count_N')
                            )
                        ),
                        html.H6("Junctions"),
                        dash_table.DataTable(
                            id=app_prefix + "junction_table",
                            columns=[
                                dict(id="sample", name="sample", type="text"),
                                dict(id="gene", name="gene", type="text"),
                                dict(id="chr", name="chr", type="numeric"),
                                dict(id="start", name="start", type="numeric"),
                                dict(id="end", name="end", type="numeric"),
                                dict(id="split_read_count", name="split_read_count", type="numeric"),
                                dict(id="split_pair_count", name="split_pair_count", type="numeric"),
                                dict(id="ds_coverage", name="ds_coverage", type="text"),
                                dict(id="us_coverage", name="us_coverage", type="text"),
                            ],
                            #data=T_subset.to_dict('records'),  #TODO:
                            page_size=10,
                            fixed_rows={'headers':True},
                            editable=False,
                            cell_selectable=True,
                            sort_action="native",
                            filter_action="native",
                            #export_columns="all",
                            export_format="csv",
                            export_headers="names",
                            style_data={
                                'width': '130px', 'minWidth': '130px', 'maxWidth': '130px',
                                'overflow': 'hidden',
                                'textOverflow': 'ellipsis',
                            },
                        ),
                    ],
                    id="right-column",
                    className="pretty_container nine columns",
                ),
            ],
            className="row flex-display",
        ),
        html.Div(
            [
                html.Div(
                    [
                        dcc.Graph(
                            id=app_prefix + "split_read_count_summary",
                            figure=px.scatter(
                                get_sample_summary_table(tumor_type),
                                x="split_read_count_N",
                                y="split_read_count_T",
                                title="Tumor split read count vs. Normal split read count of samples",
                                hover_data=["sample","gene_count_N","gene_count_T"],
                                labels={
                                    "split_read_count_N": "split read count (normal)",
                                    "split_read_count_T": "split read count (tumor)",
                                    "gene_count_N": "gene count (normal)",
                                    "gene_count_T": "gene count (tumor)",
                                },
                                trendline="ols",
                                #trendline_color_override="",
                            )
                        ),
                        dcc.Graph(
                            id=app_prefix + "split_read_count_summary_histogram",
                            figure=ff.create_distplot(
                                [
                                    get_sample_summary_table(tumor_type)["split_read_count_T"],
                                    get_sample_summary_table(tumor_type)["split_read_count_N"]
                                ],
                                ["Tumor","Normal"],
                                bin_size=20,
                                histnorm="probability",
                                curve_type="kde",
                            ).update_layout(title="split read distribution")
                        ),
                        dcc.Graph(
                            id=app_prefix + "gene_count_summary",
                            figure=px.scatter(
                                get_sample_summary_table(tumor_type),
                                x="gene_count_N",
                                y="gene_count_T",
                                title="Tumor gene count vs. Normal gene count of samples",
                                hover_data=["sample","split_read_count_N","split_read_count_T"],
                                labels={
                                    "split_read_count_N": "split read count (normal)",
                                    "split_read_count_T": "split read count (tumor)",
                                    "gene_count_N": "gene count (normal)",
                                    "gene_count_T": "gene count (tumor)",
                                },
                                trendline="ols",
                                #trendline_color_override="",
                            )
                        ),
                        dcc.Graph(
                            id=app_prefix + "gene_count_summary_histogram",
                            figure=ff.create_distplot(
                                [
                                    get_sample_summary_table(tumor_type)["gene_count_T"],
                                    get_sample_summary_table(tumor_type)["gene_count_N"]
                                ],
                                ["Tumor","Normal"],
                                bin_size=10,
                                histnorm="probability",
                                curve_type="kde",
                            ).update_layout(title="gene distribution")
                        ),
                    ],
                    className="pretty_container twelve columns",
                ),
            ],
            className="row flex-display",
        ),
        html.Div(
            [
                html.Div(
                    [
                        dcc.Graph(
                            id=app_prefix + "sample_count_summary",
                            figure=px.scatter(
                                get_gene_summary_table(tumor_type),
                                x="sample_count_N",
                                y="sample_count_T",
                                title="Tumor sample count vs. Normal sample count of genes",
                                hover_data=["gene","split_read_count_N","split_read_count_T"],
                                labels={
                                    "split_read_count_N": "split read count (normal)",
                                    "split_read_count_T": "split read count (tumor)",
                                    "gene_count_N": "gene count (normal)",
                                    "gene_count_T": "gene count (tumor)",
                                },
                                trendline="ols",
                                #trendline_color_override="",
                            )
                        ),
                        dcc.Graph(
                            id=app_prefix + "sample_count_summary_histogram",
                            figure=ff.create_distplot(
                                [
                                    get_gene_summary_table(tumor_type)["sample_count_T"],
                                    get_gene_summary_table(tumor_type)["sample_count_N"]
                                ],
                                ["Tumor","Normal"],
                                bin_size=10,
                                histnorm="probability",
                                curve_type="kde",
                            ).update_layout(title="sample count distribution")
                        ),
                    ],
                    className="pretty_container twelve columns",
                ),
            ],
            className="row flex-display",
        ),
    ],
    id="mainContainer",
    style={"display": "flex", "flex-direction": "column"},
)

if __name__ == '__main__':
    app.run_server(debug=True)
