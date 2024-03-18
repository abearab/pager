import os 
import pandas as pd 


def merge_onePAGE_results(path_list,merge_by=None,fix_bins=False,clean_geneset_names=None):
    
    genesets = []
    matrix_list = []
    summary_list = []
    
    for path in path_list:
        exp = os.path.basename(os.path.normpath(path)).split('_onePAGE_')[0]
        geneset = os.path.basename(os.path.normpath(path)).split('_onePAGE_')[1]
            
        if merge_by=='exp':
            genesets.append(exp)
        else:
            genesets.append(geneset)
        # read_onePAGE_matrices
        matrix_list.append(pd.read_csv(f'{path}/{exp}.txt.matrix',sep='\t'))
        summary_list.append(pd.read_csv(f'{path}/{exp}.txt.summary',sep='\t'))
    
    if clean_geneset_names:
        for pat in clean_geneset_names:
            genesets = [geneset.replace(pat,'') for geneset in genesets]
    
    if fix_bins:
        new_columns = matrix_list[0].columns
        matrix_list = [
            m.rename(columns=dict([
                (co,cp) for co,cp in zip(m.columns,new_columns)
            ]))
            for m in matrix_list
        ]
        matrix  = pd.concat(matrix_list)
        matrix['MOTIF'] = genesets
    
    else:
        matrix  = pd.concat(matrix_list)
        matrix['MOTIF'] = genesets
    
    summary = pd.concat(summary_list)
    summary['index'] = genesets

    return matrix,summary


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
