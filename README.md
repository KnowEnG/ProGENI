# *ProGENI* - Prioritization of Genes Enhanced with Network Information
#### Amin Emad
#### KnowEnG BD2K Center of Excellence
#### University of Illinois Urbana-Champaign

# Motivation
Identification of genes whose basal mRNA expression can predict a phenotype (e.g. sensitivity of tumor cells to different treatments) in an important task in bioinformatics. For example, in the context of genes predictors of drug sensitivity, screening the expression of these genes in the tumor tissue may suggest the best course of chemotherapy or suggest a combination of drugs to overcome chemoresistance. We present ProGENI (Prioritization of Genes Enhanced with Network Information) a novel computational method to identify such genes by leveraging their basal expressions and prior knowledge in the form of an interaction network. ProGENI is based on identifying a small set of genes where a combination of their expression and the activity level of the network module surrounding them shows a high correlation with drug response, followed by the ranking of the genes based on their relevance to this set using random walk techniques.

The figure below depcits the method overview in the context of drug response. 

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
- Gene expression (features) file: This is a genes x samples csv file where the first column contains name of genes and the first row contains name/IDs of the samples. ProGENI assumes that the expression of each gene (across all samples) follows a normal distribution. As a result, we recommend you perform proper transformation on your expression data (e.g. log2 transform on microarray data) to satsify this condition for best results. NAs are not allowed in this file. 
- Phenotype (response) file: This is a phenotype x samples csv file where the first column contains name of different phenotypes (e.g. different drugs) and the first row contains name/IDs of the samples. Make sure that the samples are ordered the same as the gene expression file (i.e. the first row of both files hould be identical). NAs are allowed in this file. 
- Network file: This is a csv file which contains information on gene-gene interactions. The first should be the header of the file. The network should be represented as a three-column format where each edge in the network is represented as a row in the file: the first two columns contain name of genes and the third column shows the weight (e.g. representing confidence) corresponding to this relationship. For example an edge between gene_1 and gene_2 can be represented as a row 'gene_1, gene_2, 0.6' in the file. If your network is not weighted, simply put 1 in the third column for all edges. 
