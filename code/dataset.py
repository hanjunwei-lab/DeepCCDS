import os
import random
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset



class NPweightingDataSet(Dataset):
    def __init__(self,
             response_fpath, 
             drug_input, 
             exprs_input,
             mut_input,
             response_type='IC50'):
        
        self.response_df = pd.read_csv(response_fpath, sep=',', header=0)
        self.drug_fp_df = pd.read_csv(drug_input, sep=',', header=0)
        self.cell_exprs_df = pd.read_csv(exprs_input, sep=',', header=0)
        self.mut_score_df = pd.read_csv(mut_input, sep=',', header=0)
        self.response_type = response_type
    
    def __len__(self):
        return self.response_df.shape[0]
    
    def __getitem__(self, index):
        response_sample = self.response_df.iloc[index]
        cell_idx = response_sample['cell_idx']
        drug_idx = response_sample['drug_idx']
        
        response_val = torch.tensor(response_sample[self.response_type]).float()
        cell_vec = torch.from_numpy(np.array(self.cell_exprs_df[self.cell_exprs_df['cell_idx'] == cell_idx].iloc[:,2:])).float().squeeze(0)
        drug_vec = torch.from_numpy(np.array(self.drug_fp_df[self.drug_fp_df['drug_idx'] == drug_idx].iloc[:,2:])).float().squeeze(0)
        mut_score = torch.from_numpy(np.array(self.mut_score_df[self.mut_score_df['cell_idx'] == cell_idx].iloc[:,2:])).float().squeeze(0)
        
        return (drug_vec,cell_vec,mut_score), response_val

