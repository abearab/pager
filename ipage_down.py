import os
import re
import gzip
from glob import glob
import itertools
import pandas as pd


def read_page_index(gs, PATH):
    """
    Find genes associated to given geneset in *index.txt files 
    """

    def read_it(raw):
        lines = [line for line in raw.read().splitlines()]
        lines = [line.split('\t') for line in lines if re.search(gs, line)]

        genes = [line[0] for line in lines]

        return genes

    if PATH.endswith('.gz'):
        with gzip.open(PATH, 'rt') as raw:
            genes = read_it(raw)

    else:
        with open(PATH) as raw:
            genes = read_it(raw)

    return genes


def read_page_names(gs, PATH):
    """
    Find annotations associated to given geneset in *names.txt files 
    """

    def read_it(raw):
        lines = [line for line in raw.read().splitlines()]
        anns = [line.split('\t') for line in lines if re.search(gs, line)]

        return anns

    if PATH.endswith('.gz'):
        with gzip.open(PATH, 'rt') as raw:
            anns = read_it(raw)

    else:
        with open(PATH) as raw:
            anns = read_it(raw)

    return anns


def read_page_annotations(gs, gs_clst, ANNDIR, gz):
    """
    Read geneset annotations and list of genes from PAGE_DATA format (*index.txt and *names.txt files)
    """
    annotations = {}

    if gz:
        names_path = glob(f'{ANNDIR}/{gs_clst}/*_names.txt.gz')
        index_path = glob(f'{ANNDIR}/{gs_clst}/*_index.txt.gz')
    else:
        names_path = glob(f'{ANNDIR}/{gs_clst}/*_names.txt')
        index_path = glob(f'{ANNDIR}/{gs_clst}/*_index.txt')

    genes = read_page_index(gs, index_path[0])
    anns = read_page_names(gs, names_path[0])

    annotations[anns[0][0]] = {}
    annotations[anns[0][0]]['names'] = anns[0][1:] + [gs_clst]
    annotations[anns[0][0]]['genes'] = genes

    return annotations


def read_pvmatrix(PATH):
    """
    Read pvmatrix.txt files into dataframes
    """
    pv_df = pd.read_csv(PATH, sep='\t', index_col=0)
    genesets = [geneset.split(' ')[0] for geneset in pv_df.index.tolist()]
    pv_df.index = genesets

    return pv_df


def style_clean_pvmatrix(pv_df):
    a, b = pv_df.columns
    return pv_df.sort_values([a, b], ascending=[True, False]).style.applymap(
        lambda x: "background-color: red" if x > 2 else "background-color: white"
    )


def read_pvmatrix_killed(PATH):
    """
    Read pvmatrix.txt.killed files that contains list of related gene sets those enriched (pvmatrix.txt)
    """
    with open(PATH) as killed:
        lines = [line for line in killed.read().splitlines()]
        out = {}
        for line in lines:
            if '\t' not in line:
                geneset0 = line.split(', ')[0]
                out[geneset0] = []
            else:
                geneset1 = line.split(', ')[0].replace('\t', '')
                out[geneset0].append(geneset1)

    return out


def make_ipage_run_data_frame(parent_dir, clean=True):
    """
    some note
    """
    gs_clusters = [os.path.basename(os.path.normpath(d)) for d in glob(parent_dir + '*/')]
    pv_dfs = []
    cols = read_pvmatrix(f'{parent_dir}/{gs_clusters[0]}/pvmatrix.txt').columns.to_list()
    for gs_cluster in gs_clusters:
        pv_df = read_pvmatrix(f'{parent_dir}/{gs_cluster}/pvmatrix.txt')
        pv_df.insert(0, 'gs_cluster', gs_cluster)
        pv_df = pv_df.set_index([pv_df.index, pv_df.gs_cluster]).rename_axis(['gene_set', 'gs_cluster'])
        pv_df = pv_df.drop(columns=['gs_cluster'])
        pv_df.columns = cols
        if not clean:
            pv_dfs.append(pv_df)
        else:
            pv_dfs.append(pv_df[(pv_df.iloc[:, 0] >= 2) | (pv_df.iloc[:, 10] >= 2)].iloc[:, [0, 10]])

    out = pd.concat(pv_dfs, axis=0)
    return out


def make_annotation_dict(pv_df, ANNDIR=f'{os.path.dirname(__file__)}/annotations/', gz=True, verbose=False):
    """
    Read annotations for given pv dataframe(s)
    """
    ann = {}

    if type(pv_df) is pd.core.frame.DataFrame:
        sets = set(sorted(pv_df.index.to_list(), key=lambda x: x[0]))

    elif type(pv_df) is list:
        raw = []
        for df in pv_df:
            raw = raw + df.index.to_list()
        sets = set(sorted(raw, key=lambda x: x[0]))

    for gene_set, gs_cluster in sets:
        if verbose:
            print(gs_cluster, gene_set)
        ann.update(read_page_annotations(gene_set, gs_cluster, ANNDIR, gz))

    return ann


# def check_gene(ann,genes):
#     membs = {}
#     for gs in ann:
#         membs[gs]= []
#         for gene in genes:
#             if gene in ann[gs]['genes']:
#                 membs[gs].append(gene)
#     return membs


# def make_ipage_run_dict(parent_dir, ANNDIR=f'{os.path.dirname(__file__)}/annotations/', gz=True):
#     """
#     Include annotations for the genesets in each result tables
#     """
#     out = {}
#
#     pv_df = make_ipage_run_data_frame(parent_dir)
#     out['pv_df'] = pv_df
#     # read_pvmatrix_killed
#
#     ann = make_annotation_dict(pv_df, ANNDIR, gz)
#     out['ann'] = ann
#     return out
