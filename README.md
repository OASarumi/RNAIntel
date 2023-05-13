# RNAIntels
Machine learning-based tool for discovering coding and non coding small RNAs in human RNA-seq data.
## Description
Following receiving raw sequence readings from an organism, next is the need to classify them into their functional units–coding or non-coding sequences. Coding RNA sequences are translated into protein, while non-coding RNAs, e.g., small non-coding RNA such as miRNA or siRNA, are not translated into proteins but are essential for many significant biological processes, such as gene regulation and gene expressions. Although several research efforts are geared towards developing computational tools for functional prediction of RNA-seq data, identifying coding and non-coding RNAs from raw RNA-seq data remains a challenging task. Here we present a software tool using machine-learning approaches to classify short, human RNA sequences (in the range of 28 – 51 nucleotides) into coding and non-coding RNAs. In contrast to many existing approaches, RNAIntels neither requires raw RNA-Seq sequences (FASTQ) nor aligning the sequences to a human reference genome, thus reducing complexities issues that made other tools less attractive to users.
## Requirements
#### Hardware Requirements
RNAIntels requires a  functional computer system with a minimum of 4g RAM for efficient encoding of the RNA-seq data.

#### OS Requirements
RNAIntels meet the software portability requirements. The software is not platform dependent and can run on Windows, Linux and Unix operating systems. The software has been tested on Windows 10 and Unix operating systems


## Dependencies
The code is written in Python 3.9.16 and tested on Unix OS with the following libraries installed:

Library | Version
--- | --- 
tensorflow | 2.10.0
biopython | 1.78
pandas | 1.53
matplotlib | 3.7.1
numpy | 1.23.5
scikit-learn| 1.2.2
keras | 2.10.0
plotly | 5.9.0


## Data
The training data, GRCh38 homo sapiens reference data for building the machine learning model, was retrieved from the ENSEMBL database (www.ensembl.org). The raw dataset contained 207,877 instances of protein-coding sequences and 63,865 instances of non-protein coding sequences. We observed that the raw data was inundated with instances of pseudogenes and overlapping sequences. After cleaning the raw data, we obtained a total of 69,420 protein-coding sequences and 28,225 non-coding sequences. We created a balanced dataset ( lenght <= 60 ) from the cleaned data  as training data for the machine learning operations. The balanced dataset can be located at Data/train_test data.

The validation data, human RNA-seq was received from the Institute for Lung Research, Universities of Giessen. This can be made available on request.

## Execution
You can execute RNAIntels by running the following codes in python enviroment. The executable file is saved as h5 format, and can be downloaded from /code/RNAIntels.h5

`from tensorflow.keras.models import load_model`

`Model = load_model ('path/to/RNAIntel.h5')`

`Model.summary()`

`y_predict =Model.predict(X_test)`
