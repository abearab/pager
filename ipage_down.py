import os
import re
import gzip
from glob import glob
import itertools
import pandas as pd


def write_page_index(input, PATH):
    """Genesets of each gene
    Write correctly formatted {PATH}_index.txt file
    input: dictionary; keys are genes and values are gene sets
    """
    with gzip.open(PATH, 'wt') as raw:
        for gn,gs in input.items():
            raw.write('\t'.join([gn] + list(gs))+'\n')


def write_page_names(input, PATH):
    """Genesets annotation
    Write {PATH}_names.txt file
    input: pandas dataframe with three columns
        1) Geneset index matched with page_index values
        2) Geneset name (human readable)
        3) Name Space from GO database (???)
    """
    input.to_csv(PATH, index=True,header=False, sep='\t', compression='gzip')
    
    
def search_page_index(gs, PATH):
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


def search_page_names(gs, PATH):
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

    genes = search_page_index(gs, index_path[0])
    anns  = search_page_names(gs, names_path[0])

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


def clean_bins_range(pv_df):
    bins = pv_df.columns.tolist()
    ranges = [bn.replace('[','').replace(']','').split(' ') for bn in bins]
    ave = sum([(float(b) - float(a)) for a,b in ranges[1:-1]]) / len(ranges[1:-1])
    ranges[0][0]  = str(round(float(ranges[0][1])  - ave,2))
    ranges[-1][1] = str(round(float(ranges[-1][0]) + ave,2))

    pv_df.columns = [f'[{a} {b}]' for a,b in ranges]
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


def read_ipage_intersections_file(gs_cluster_path,clust,gs=None):
    '''read data from `output.ipage_intersections` file from an iPAGE run
    `output.ipage_intersections` is a three column table that identify genes enriched in each cluster bins. 
    # Group	Cluster	Identifiers

    Args:
        gs_cluster_path: path to the parent directory of `output.ipage_intersections` file.
        clust: the cluster (bin) number to extract identifiers.
        gs: gene set name to search .
    Returns:
        a dictionary in which keys are pathway names and values are genes reported in given `clust`.
    '''
    # Group	Cluster	Identifiers
    with open(f'{gs_cluster_path}/output.ipage_intersections') as raw:
        lines = [line for line in raw.read().splitlines()]
        if gs: 
            lines = [line.split('\t') for line in lines if re.search(gs, line)]
        else:
            lines = [line.split('\t') for line in lines]

    return dict([(line[0].split(' ')[0],line[2:]) for line in lines if line[1] == clust ])


def merge_multiple_pvmat(pvmat_list):
    '''merge multiple `pvmatrix.txt` files from different iPAGE runs (on one input) 

    Args:
        pvmat_list: list of pvmatrix file path.
    Returns:
        pandas data.frame
    '''
    df = ipd.clean_bins_range(
        ipd.read_pvmatrix(pvmat_list[0])
    )
    
    cols = df.columns
    
    df = pd.concat(
        [df] + [
            ipd.read_pvmatrix(pvmat).set_axis(cols, axis=1, inplace=False) 
            for pvmat in pvmat_list[1:]
        ]
    )
    
    df = df.groupby(df.index).first()
    # ipd.style_clean_pvmatrix(df.iloc[:,[0,10]])    
    return df


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
