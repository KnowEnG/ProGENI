# *ProGENI* - Prioritization of Genes Enhanced with Network Information
#### Amin Emad
#### KnowEnG BD2K Center of Excellence
#### University of Illinois Urbana-Champaign

# Motivation
Identification of genes whose basal mRNA expression can predict a phenotype (e.g. sensitivity of tumor cells to different treatments) in an important task in bioinformatics. For example, in the context of genes predictors of drug sensitivity, screening the expression of these genes in the tumor tissue may suggest the best course of chemotherapy or suggest a combination of drugs to overcome chemoresistance. We present ProGENI (Prioritization of Genes Enhanced with Network Information) a novel computational method to identify such genes by leveraging their basal expressions and prior knowledge in the form of an interaction network. ProGENI is based on identifying a small set of genes where a combination of their expression and the activity level of the network module surrounding them shows a high correlation with drug response, followed by the ranking of the genes based on their relevance to this set using random walk techniques.

The figure below depcits the method overview in the context of drug response. (a) shows an overview of the ProGENI method, while (b) shows how bootstrap sampling can be used in order to obtain robust ranking. 

![Method Overview](/images/Pipeline_ProGENI.png)

# Requirements

In order to run the code, you need to have Python 3.5 installed. In addition, the code uses the following python modules/libraries which need to be installed:
- [Numpy](http://www.numpy.org/)
- [Scipy](https://www.scipy.org/)
- [Pandas](http://pandas.pydata.org/)
- [Sklearn (scikit-learn)](http://scikit-learn.org/stable/)

Instead of installing all these libraries independently, you can use prebulit Python distributions such as [Anaconda](https://www.continuum.io/downloads), which provides a free academic subscription.

# Input files

ProGENI requires three files as input:
#### Gene expression (features) file:
This is a genes x samples csv file where the first column contains name of genes and the first row contains name/IDs of the samples. ProGENI assumes that the expression of each gene (across all samples) follows a normal distribution. As a result, we recommend you perform proper transformation on your expression data (e.g. log2 transform on microarray data) to satsify this condition for best results. NAs are not allowed in this file. 

Example Gene expression file:

|  | sample_1 | sample_2 | sample_3 |
| :--- | :--- | :--- | :--- |
| G1 | 0.24 | 0.67 | 2.12 |  
| G2 | 0.34 | -1.34 | 0.45 |
| G3 | 1.51 | 0.05 | -0.22 |
| G4 | 0.03 | 0.55 | 1.15 |
| G5 | -0.23 | 0.23 | 0.55 |
| G6 | 0.94 | 0.33 | 1.12 |


#### Phenotype (response) file:
This is a phenotype x samples csv file where the first column contains name of different phenotypes (e.g. different drugs) and the first row contains name/IDs of the samples. Make sure that the samples are ordered the same as the gene expression file (i.e. the first row of both files hould be identical). NAs are allowed in this file. 

#### Network edge file:
This is a csv file which contains information on gene-gene interactions. The first should be the header of the file. The network should be represented as a three-column format where each edge in the network is represented as a row in the file: the first two columns contain name of genes and the third column shows the weight (e.g. representing confidence) corresponding to this relationship. If the set of genes in the network is slightly different from the set of genes in the gene expression data, ProGENI will focus on the intersection of the genes.  

Example network edge file:
```
node_1  node_2  weight
G1	G4	0.76
G1      G6      0.24
G2	G1	0.45
G2	G7	0.55
G4	G5	0.23
```

