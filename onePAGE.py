import os 
import pandas as pd 

def merge_onePAGE_results(path_list):
    path0 = path_list[0]
    
    pdir = '/'.join(path0.split('/')[:-1])
    exp = path0.split('/')[-1].split('_onePAGE_')[0]
    onePAGE = f'{pdir}/{exp}_onePAGE'
    if not os.path.exists(onePAGE): 
        os.mkdir(onePAGE)
    if not os.path.exists(f'{pdir}/{exp}_onePAGE/{exp}.txt'):
        os.rename(f'{path0}/{exp}.txt', f'{pdir}/{exp}_onePAGE/{exp}.txt')
    
    genesets = []
    matrix_list = []
    summary_list = []
    
    for path in path_list:
        geneset = path.split('/')[-1].split('_onePAGE_')[1]
        genesets.append(geneset)
        # read_onePAGE_matrices
        matrix_list.append(pd.read_csv(f'{path}/{exp}.txt.matrix',sep='\t'))
        summary_list.append(pd.read_csv(f'{path}/{exp}.txt.summary',sep='\t'))
    
    summary = pd.concat(summary_list)
    summary['index'] = genesets
    matrix  = pd.concat(matrix_list)
    matrix['MOTIF'] = genesets

    matrix.to_csv(f'{onePAGE}/{exp}.txt.matrix',sep='\t',header=True,index=False)
    summary.to_csv(f'{onePAGE}/{exp}.txt.summary',sep='\t',header=True,index=False)
    
    return matrix


def merge_onePAGE_results_1(path_list,return_df=False):
    path0 = path_list[0]
    
    pdir = '/'.join(path0.split('/')[:-1])
    geneset = path0.split('/')[-1].split('_onePAGE_')[1]
    merge_id = path0.split('/')[-1].split('_onePAGE_')[1]
    onePAGE = f'{pdir}/{geneset}_onePAGE'
    if not os.path.exists(onePAGE): 
        os.mkdir(onePAGE)
    # if not os.path.exists(f'{pdir}/{exp}_onePAGE/{exp}.txt'):
    #     os.rename(f'{path0}/{exp}.txt', f'{pdir}/{merge_id}/{exp}.txt')
    
    samples = []
    matrix_list = []
    summary_list = []
    
    for path in path_list:
        exp = path.split('/')[-1].split('_onePAGE_')[0]
        samples.append(exp)
        # read_onePAGE_matrices
        matrix_list.append(pd.read_csv(f'{path}/{exp}.txt.matrix',sep='\t'))
        summary_list.append(pd.read_csv(f'{path}/{exp}.txt.summary',sep='\t'))
    
    summary = pd.concat(summary_list)
    summary['index'] = samples
    matrix  = pd.concat(matrix_list)
    matrix['MOTIF'] = samples

    if return_df:
        return matrix, summary
    else:
        matrix.to_csv(f'{onePAGE}/{geneset}.txt.matrix',sep='\t',header=True,index=False)
        summary.to_csv(f'{onePAGE}/{geneset}.txt.summary',sep='\t',header=True,index=False)
