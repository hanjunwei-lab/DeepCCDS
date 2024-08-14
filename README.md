# DeepCCDS - Interpretable deep learning framework for predicting cancer cell drug sensitivity through characterizing cancer driver signals
-----------------------------------------------------------------
Code by **Jiashuo Wu** and **Jiyin Lai** at Harbin Medical University.

## 1. Introduction
This repository contains source code and data for **DeepCCDS** 

**DeepCCDS** is a **Deep** learning framework for **C**ancer **C**ell **D**rug **S**ensitivity prediction through **C**haracterizing **C**ancer **D**river **S**ignals. **DeepCCDS** incorporates a prior knowledge network to characterize cancer driver signals, building upon the self-supervised neural network framework. The signals can reflect key mechanisms influencing cancer cell development and drug response, enhancing the model's predictive performance and interpretability.

**DeepCCDS** demonstrated superior predictive performance compared to other state-of-the-art methods across multiple datasets. Benefiting from integrating prior knowledge, DeepCCDS exhibits powerful feature representation capabilities and interpretability. Based on these feature representations, we have identified embedding features that could potentially be used for drug screening in new indications. Further, we demonstrate the applicability of DeepCCDS on solid tumor samples from The Cancer Genome Atlas. We believe integrating DeepCCDS into clinical decision-making processes can potentially improve the selection of personalized treatment strategies for cancer patients.

## 2. Design of DeepCCDS

![alt text](../image/Overall_architecture.tif "DeepCCDS")

Figure 1: Overall architecture of DeepCCDS

## Installation

**DeepCCDS** relies on [R (version 4.3.1)](https://cran.r-project.org/bin/windows/base/old/4.3.1/) and [Python (version 3.10.12)](https://www.python.org/downloads/release/python-31012/) environments. Before using DeepCCDS, you must set up these two environments and install some necessary packages.

- Install the necessary R packages for **DeepCCDS**:

```sh
install.packages(c("dplyr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GSVA", "clusterProfiler", "org.Hs.eg.db"))
```
- Install the necessary Python packages for **DeepCCDS**:

```sh
conda install numpy pandas scipy scikit-learn pytorch
```

## 3. Before usage 
We have made available all the code necessary to execute **DeepCCDS** in this GitHub repository. Please ensure that you replace the input data paths in the code with your own storage locations.

## 4. Usage

The implemention of **DeepCCDS** framework is divided into two main parts: 

- **Part I**: Characterizing cancer driver signals as pathway representations through a prior knowledge network.
- **Part II**: The complete training stage, in which the entire DeepCCDS framework undergoes end-to-end training.

### 4.1. Implemention of Part I
**1. Code**: The code used to implement **Part I** is in the file "[Characterizing_CDS.R](code/Characterizing_CDS.R)" which located at folder ``code/``. This file encompasses the original data processing, execution of the random walk algorithm and enrichment analysis.

**2. Data**: The datasets used to implement **Part I** can be downloaded from the ``scDEAL.7z`` link:
[Click here to download scDEAL.zip from google drive](https://drive.google.com/file/d/1egI0B5YiDrHQz-4jiWClfinQqVqNbvF-/view?usp=drive_link)

The file ``scDEAL.7z`` includes all the input datasets of ``Characterizing_CDS.R``:

| File                              | Description                                                                   |
|------------------------------------|------------------------------------------------------------------------|
| [model_list_20230923](https://cellmodelpassports.sanger.ac.uk/downloads)                             | Cell line annotation information of GDSC                            |
| [mutations_all_20230202](https://cellmodelpassports.sanger.ac.uk/downloads)                            | Cell line mutation information of GDSC |
| [rnaseq_tpm_20220624](https://cellmodelpassports.sanger.ac.uk/downloads)                           | Cell line transcription information of GDSC                              |
| [gene_identifiers_20191101](https://cellmodelpassports.sanger.ac.uk/downloads) | Correspondence of different identifiers of genes                                       |
| [Census_allMon](https://cancer.sanger.ac.uk/census)                 | Cancer drivers from CGC                              |
| [adjM_PPI](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01119-6
        
        )                           | PPI networks obtained from previous studies          |
| [listp243](https://www.genome.jp/kegg/)                           | Kegg pathways and excluded all disease pathways                |

**3. Output**: The final output of the file ``Characterizing_CDS.R`` is a pathway activity matrix that is used for subsequent training of the complete model.

### 4.2. Implemention of Part II  

**1. Code**: The complete training code for **DeepCCDS** is in the file ``main_DeepCCDS.ipynb`` which located at folder ``code/``. The file include some tunable parameters.  

   - ``--learning_rate``:       learning_rate value
   - ``--epochs``:         number of training iterations
   - ``--batch_size``:        number of samples per batch
   - ``--code_dropout_rate``:  dropout rate value
   - ``--weight_decay``:   L2 regularization coefficient        
   - ``--patience``: patience of the early stop
   - ``--testset_yes``:   whether to use test sets    
   - ``--code_dim``:     self-encoder bottleneck layer dimension
   - ``--drug_hidden_dims``:    drug self-encoder hidden layer dimension
   - ``--mut_hidden_dims``:     mutation self-encoder hidden layer dimension
   - ``--forward_net_hidden_dim1``:     feedforward neural network hidden layer dimension
   - ``--forward_net_hidden_dim2``:     feedforward neural network hidden layer dimension   


**2. Data**: The datasets used to train **DeepCCDS** are located at folder ``data/``
| File                              | Description                                                                   |
|------------------------------------|------------------------------------------------------------------------|
| GDSC_mutation_input.csv                             | Mutation features of cell lines                            |
| GDSC_ssgsea_input.csv                           | Pathway activities of cell lines|
| GDSC_SMILE_input.csv                           | Drug structure features of cell lines                              |
| GDSC_train_IC50_by_borh_cv00.csv | Training samples (cell-drug pairs)                                       |
| GDSC_train_IC50_by_borh_cv00.csv                 | Validating samples (cell-drug pairs)                              |
| GDSC_test_IC50_by_borh_cv00.csv                           | Testing samples (cell-drug pairs)          |

**3. Output**: Get the best trained model. The models used in this paper are located at folder ``model/``.

## 5. Feature attribution analysis

To explore the relationship between each cell embedded features and drug sensitivity, we employed the **Integrated Gradients (IG)** method ([https://doi.org/10.48550/arXiv.1703.01365](https://doi.org/10.48550/arXiv.1703.01365)). **IG** attributes the model's prediction for its input features by computing gradients for each input and measures the change in the output based on the small changes in the input. We calculated the average attribution of features across all samples to represent the global importance, termed the **IG score**. The calculation was performed through the ``IntegratedGradients`` class from the Python ``Captum`` library.

Install the package:

```sh
conda install Captum
```
The code for calculating **IG score** is in the file ``feature_attribution.ipynb`` which located at folder ``code/``.
