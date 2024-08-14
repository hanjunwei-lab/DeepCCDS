# === Fixed === #
import os
import random
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")
from collections import defaultdict
import torch
import torch.nn as nn
import torch.optim as optim
from torch.optim import lr_scheduler

from sklearn.metrics import mean_squared_error
from lifelines.utils import concordance_index
from scipy.stats import pearsonr,spearmanr
import scipy.stats as stats

# === Task-specific === #
import sys
sys.path.append('..')
from copy import deepcopy
from dataset import NPweightingDataSet
from utils import *
from trainers import logging, train, validate, test

# ====== Random Seed Initialization ====== #
def seed_everything(seed = 3078):
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = True
seed_everything()

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)


# ====== Model Definition ====== #
class DeepAutoencoderThreeHiddenLayers(nn.Module):
    def __init__(self, input_dim, hidden_dims, code_dim, activation_func=nn.ReLU,
                 code_activation=True, dropout=False, dropout_rate=0.5):
        super(DeepAutoencoderThreeHiddenLayers, self).__init__()
        # Establish encoder
        modules = []
        modules.append(nn.Linear(input_dim, hidden_dims[0]))
        modules.append(activation_func())
        if dropout:
            modules.append(nn.Dropout(dropout_rate))

        
            
        for input_size, output_size in zip(hidden_dims, hidden_dims[1:]):
            
            modules.append(nn.Linear(input_size, output_size))
            modules.append(activation_func())
            if dropout:
                modules.append(nn.Dropout(dropout_rate))
                
        modules.append(nn.Linear(hidden_dims[-1], code_dim))
        if code_activation:
            modules.append(activation_func())
        self.encoder = nn.Sequential(*modules)
        
        # Establish decoder
        modules = []
        
        modules.append(nn.Linear(code_dim, hidden_dims[-1]))
        modules.append(activation_func())
        if dropout:
            modules.append(nn.Dropout(dropout_rate))
        
        for input_size, output_size in zip(hidden_dims[::-1], hidden_dims[-2::-1]):
            modules.append(nn.Linear(input_size, output_size))
            modules.append(activation_func())
            if dropout:
                modules.append(nn.Dropout(dropout_rate))
        modules.append(nn.Linear(hidden_dims[0], input_dim))
        modules.append(nn.Sigmoid())
        self.decoder = nn.Sequential(*modules)
        
    def forward(self, x):
        x = self.encoder(x)
        code = x
        x = self.decoder(x)
        return code, x
    

class ForwardNetworkTwoHiddenLayers(torch.nn.Module):
    def __init__(self, input_dim, hidden_dim1, hidden_dim2, activation_func=nn.ReLU,
                out_activation=None):
        super(ForwardNetworkTwoHiddenLayers, self).__init__()
        
        self.layers = nn.Sequential(
             nn.Linear(input_dim, hidden_dim1),
             nn.BatchNorm1d(hidden_dim1),
             activation_func(),
             nn.Linear(hidden_dim1, hidden_dim2),
             nn.BatchNorm1d(hidden_dim2),
             activation_func(),
             nn.Linear(hidden_dim2, 1))
        
        self.out_activation = out_activation
        
    
    def forward(self, x):
        if self.out_activation:
            return self.out_activation(self.layers(x))
        else:
            return self.layers(x)

class DEERS_Concat(torch.nn.Module):
    def __init__(self, drug_autoencoder, mut_line_autoencoder,
                 forward_network):
        super(DEERS_Concat, self).__init__()
        self.drug_autoencoder = drug_autoencoder
        self.mut_line_autoencoder = mut_line_autoencoder
        self.forward_network = forward_network
        
    def forward(self, drug_features, mut_features, cell_features):
        
        drug_code, drug_reconstruction = self.drug_autoencoder(drug_features)
        mut_code, mut_reconstruction = self.mut_line_autoencoder(mut_features)
        x = torch.cat((drug_code, mut_code, cell_features), axis=1)    ####
        return self.forward_network(x), drug_reconstruction, mut_reconstruction
    
class MergedLoss(nn.Module):
    def __init__(self, y_loss_weight=1., drug_reconstruction_loss_weight=0.1, mut_reconstruction_loss_weight=0.2):
        super(MergedLoss, self).__init__()
        self.y_loss_weight = y_loss_weight
        self.drug_reconstruction_loss_weight = drug_reconstruction_loss_weight
        self.mut_reconstruction_loss_weight = mut_reconstruction_loss_weight
        self.output_criterion = nn.MSELoss()
        self.reconstruction_criterion = nn.BCELoss()
    
    def forward(self, pred_y, drug_reconstruction, mut_reconstruction,drug_input,mut_input, true_y):
        output_loss = self.output_criterion(pred_y, true_y)
        drug_reconstruction_loss = self.reconstruction_criterion(drug_reconstruction, drug_input)
        mut_reconstruction_loss = self.reconstruction_criterion(mut_reconstruction, mut_input)
        return output_loss, drug_reconstruction_loss,mut_reconstruction_loss