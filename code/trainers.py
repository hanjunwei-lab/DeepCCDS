import time
import os
import torch
from torch import nn
import torch.nn.functional as F
import numpy as np
from sklearn.metrics import mean_squared_error
from lifelines.utils import concordance_index
from scipy.stats import pearsonr,spearmanr


    
def logging(msg, outdir, log_fpath):
    fpath = os.path.join(outdir, log_fpath)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    with open(fpath, 'a') as fw:
        fw.write("%s\n" % msg)
    print(msg)


def train(model, epoch, train_loader, optimizer, loss_fn, device):
    # ====== Train ====== #
   
    list_train_batch_loss = []
    list_train_batch_out  = []
    list_train_batch_true = []
    
    model.train()

    for batch_idx, (data, target) in enumerate(train_loader):
        print(batch_idx)
        drug_input, cell_input, mut_input = tuple(d.to(device) for d in data)
        true_y = target.unsqueeze(1).to(device)
        
        optimizer.zero_grad()
        pred_y, drug_reconstruction, mut_reconstruction= model(drug_input,mut_input,cell_input)
        output_loss, drug_reconstruction_loss, mut_reconstruction_loss = loss_fn(pred_y, drug_reconstruction, mut_reconstruction,drug_input, mut_input,true_y)
        
        loss = loss_fn.y_loss_weight * output_loss + \
               loss_fn.drug_reconstruction_loss_weight * drug_reconstruction_loss + \
               loss_fn.mut_reconstruction_loss_weight * mut_reconstruction_loss
        
        loss.backward()
        optimizer.step()
        
        list_train_batch_out.extend(pred_y.detach().cpu().numpy())
        list_train_batch_true.extend(true_y.detach().cpu().numpy())
        list_train_batch_loss.append(loss.detach().cpu().numpy())
        
    return model, list_train_batch_loss, list_train_batch_out, list_train_batch_true


def validate(model, valid_loader, loss_fn, device):
    list_val_batch_loss = []
    list_val_batch_out  = []
    list_val_batch_true = []
    
    model.eval()
    for batch_idx, (data, target) in enumerate(valid_loader):
        drug_input, cell_input, mut_input = tuple(d.to(device) for d in data)
        true_y = target.unsqueeze(1).to(device)
        pred_y, drug_reconstruction, mut_reconstruction = model(drug_input,mut_input,cell_input)
        output_loss, drug_reconstruction_loss, mut_reconstruction_loss = loss_fn(pred_y, drug_reconstruction, mut_reconstruction,drug_input, mut_input,true_y)
        
        loss = loss_fn.y_loss_weight * output_loss + \
               loss_fn.drug_reconstruction_loss_weight * drug_reconstruction_loss + \
               loss_fn.mut_reconstruction_loss_weight * mut_reconstruction_loss
        

        list_val_batch_out.extend(pred_y.detach().cpu().numpy())
        list_val_batch_true.extend(true_y.detach().cpu().numpy())

        list_val_batch_loss.append(loss.detach().cpu().numpy())
    
    return list_val_batch_loss, list_val_batch_out, list_val_batch_true





def test(model, test_loader, loss_fn, device):
    with torch.no_grad():
        # ====== Test ====== #
        list_test_loss = []
        list_test_out  = []
        list_test_true = []
        list_test_output_loss = []
        model.eval()
        for batch_idx, (data, target) in enumerate(test_loader):
            drug_input, cell_input, mut_input = tuple(d.to(device) for d in data)
            true_y = target.unsqueeze(1).to(device)
            pred_y, drug_reconstruction, mut_reconstruction = model(drug_input,mut_input,cell_input)
            output_loss, drug_reconstruction_loss, mut_reconstruction_loss= loss_fn(pred_y, drug_reconstruction, mut_reconstruction,drug_input, mut_input,true_y)
            
            loss = loss_fn.y_loss_weight * output_loss + \
                loss_fn.drug_reconstruction_loss_weight * drug_reconstruction_loss + \
                loss_fn.mut_reconstruction_loss_weight * mut_reconstruction_loss
            
            list_test_out.append(pred_y.detach().cpu().numpy().item())
            list_test_true.append(true_y.detach().cpu().numpy().item())
            list_test_loss.append(loss.detach().cpu().numpy())
            list_test_output_loss.append(output_loss.detach().cpu().numpy())
            
    test_loss= sum(list_test_loss)/len(list_test_loss)
    test_rmse = np.sqrt(mean_squared_error(np.array(list_test_out), np.array(list_test_true)))
    test_corr, _p = pearsonr(np.array(list_test_out), np.array(list_test_true))
    test_spearman, _p = spearmanr(np.array(list_test_out), np.array(list_test_true))
    test_ci = concordance_index(np.array(list_test_out), np.array(list_test_true))
    
    #print(f"Test: Loss: {test_loss:.4f}, RMSE; {test_rmse:.4f}, CORR: {test_corr:.4f}, SPEARMAN: {test_spearman:.4f}, CI: {test_ci:.4f}")
    return  test_loss,test_rmse, test_corr, test_spearman, test_ci, list_test_loss,list_test_output_loss,list_test_out, list_test_true # list_test_output_loss,


