import os 
import re 
import gzip
from glob import glob
import itertools
import pandas as pd


def read_page_index(gs,PATH):
    '''
    Find genes associated to given geneset in *index.txt files 
    '''
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
    '''
    Find annotations associated to given geneset in *names.txt files 
    '''
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


def read_page_annotations(gs,gs_clst,ANNDIR,gz):
    '''
    Read geneset annotations and list of genes from PAGE_DATA format (*index.txt and *names.txt files)
    '''
    annotations = {}
    
    if gz:
        names_path = glob(f'{ANNDIR}/{gs_clst}/*_names.txt.gz')
        index_path = glob(f'{ANNDIR}/{gs_clst}/*_index.txt.gz')
    else:
        names_path = glob(f'{ANNDIR}/{gs_clst}/*_names.txt')
        index_path = glob(f'{ANNDIR}/{gs_clst}/*_index.txt')

    genes = read_page_index(gs,index_path[0])
    anns  = read_page_names(gs,names_path[0])
    
    annotations[anns[0][0]] = {}
    annotations[anns[0][0]]['names'] = anns[0][1:]
    annotations[anns[0][0]]['genes'] = genes
    
    return annotations



def read_pvmatrix(PATH):
    '''
    Read p-value matrix (pvmatrix.txt files) into dataframes 
    '''
    pv_df = pd.read_csv(PATH, sep='\t',index_col=0) 
    genesets = [geneset.split(' ')[0] for geneset in pv_df.index.tolist()]
    pv_df.index = genesets
    
    return pv_df


# def read_pvmatrix_killed()


def make_ipage_run_dict(parent_dir, ANNDIR='/flash/bin/iPAGEv1.0/PAGE_DATA/ANNOTATIONS', gz=True):
    '''
    Include annotations for the genesets in each result tables 
    '''
    gs_clsters = [os.path.basename(os.path.normpath(d)) for d in glob(parent_dir+'*/')]
        
    out = {}
    
    # for in for would be super slow!
    for gs_clst in gs_clsters:
        pv_df = read_pvmatrix(f'{parent_dir}/{gs_clst}/pvmatrix.txt')
        ann = {}
        for gs in pv_df.index:
            ann.update(read_page_annotations(gs,gs_clst,ANNDIR,gz))
        
        out[gs_clst] = {}
        out[gs_clst]['annotations'] = ann
        out[gs_clst]['pvmatrix'] = pv_df
    
    return out
