#!/usr/bin/env python

"""
Implementation fo ProGENI
"""
import os
import sys
import argparse
import time
import warnings
import numpy as np
from numpy import mean
import pandas as pd
import scipy.stats as ss
from scipy import sparse
from scipy.stats.mstats import zscore
from scipy.stats.mstats import gmean
from sklearn.preprocessing import normalize
from scipy.sparse import SparseEfficiencyWarning

# pylint: disable=no-member
###############################################################################
def parse_args():
    """
    Parse the arguments.
    Parse the command line arguments/options using the argparse module
    and return the parsed arguments (as an argparse.Namespace object,
    as returned by argparse.parse_args()).
    Returns:
        argparse.Namespace: the parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('input_expression', type=str,
                        help='name of the file containg expression data')
    parser.add_argument('input_response', type=str,
                        help='name of the file containg response data')
    parser.add_argument('input_network', type=str,
                        help='name of the file containg network data')
    parser.add_argument('-s', '--seed', type=int, default=1011,
                        help='seed used for random generator')
    parser.add_argument('-nr', '--num_RCG', type=int, default=100,
                        help='number of genes in the response-correlated gene (RCG) set')
    parser.add_argument('-pt', '--prob_restart_trans', type=float, default=0.5,
                        help='restart probability of RWR to network-transform expression')
    parser.add_argument('-pr', '--prob_restart_rank', type=float, default=0.5,
                        help='restart probability for RWR used to rank nodes w.r.t. RCG')
    parser.add_argument('-t', '--tolerance', type=float, default=1e-8,
                        help='tolerance used to determine convergence of RWR')
    parser.add_argument('-mi', '--max_iteration', type=int, default=100,
                        help='maximum number of iterations used in RWR')
    parser.add_argument('-nb', '--num_bootstrap', type=int, default=1,
                        help='number of bootstrap samplings')
    parser.add_argument('-pb', '--percent_bootstrap', type=int, default=100,
                        help='percent of samples for bootstrap samplinga (between 0-100)')
    parser.add_argument('-de', '--directory_expression', type=str,
                        default='./',
                        help='directory containing expression data')
    parser.add_argument('-dr', '--directory_response', type=str,
                        default='./',
                        help='directory containing response data')
    parser.add_argument('-dn', '--directory_network', type=str,
                        default='./',
                        help='directory containing network data')
    parser.add_argument('-do', '--directory_out', type=str,
                        default='./',
                        help='directory for the results')
    parser.add_argument('-o', '--output', type=str,
                        default='results.csv',
                        help='name of the file containg the results')
    args = parser.parse_args()
    return args


###############################################################################
def rank_aggregate_borda(list_of_list, method):
    """
    This function receives a list of list of genes and a method. The genes are
    ranked based on Borda's method. Final output is an aggregated ranked list of 
    genes. 
    Input: 
        list_of_list: a list of list of genes
        method: "arithmetic_mean" or "geometric_mean"
    Output: 
        ranked list of genes where the first gene has the highest rank
    """
    dic_tmp = {key:[] for key in list_of_list[0]} #A dictionary with keys being gene names
    list_length = len(list_of_list[0])
    for list1 in list_of_list:
        for i in range(list_length):
            dic_tmp[list1[i]].append(list_length - i)

    if method == 'arithmetic_mean':
        dic_tmp_agg = {key:mean(dic_tmp[key]) for key in dic_tmp}
    if method == 'geometric_mean':
        dic_tmp_agg = {key:gmean(dic_tmp[key]) for key in dic_tmp}

    ranked_tmp = sorted(dic_tmp_agg, key=dic_tmp_agg.get, reverse=True)   #sort descending

    return ranked_tmp


###############################################################################
def is_number(num):
    """
    Determine whether a string s is a number (i.e., any floating point
    representation of a number, including scientific notation)
    """
    try:
        float(num)
        return True
    except ValueError:
        return False



###############################################################################
def spread_match_network(expr_df_in, node_names_in):
    """
    Matches S (spreadsheet of gene expressions) and N (network)
    The function returns expr_df_out which is formed by reshuffling columns of
    expr_df_in. Also, node_names_out is formed by reshuffling node_names_in. The
    intersection of node_names_out and column names of expr_df_out are placed at
    the beginning of both lists.
    Input: 
        expr_df_in: A pandas dataframe corresponding to gene expression 
        node_names_in: Name of the nodes in the network
    Output: 
        expr_df_out: Reorganized dataframe of gene expressions
        nodes_names_out: Reordered node names
        nodes_genes_intersect: Sorted list of shared genes
    """
    node_names_in_set = set(node_names_in)
    gene_names_in_set = set(expr_df_in.columns.values)

    nodes_genes_intersect = sorted(list(gene_names_in_set & node_names_in_set))
    nodes_minus_genes = sorted(list(node_names_in_set - gene_names_in_set))
    genes_minus_nodes = sorted(list(gene_names_in_set - node_names_in_set))

    genes_names_out = nodes_genes_intersect + genes_minus_nodes
    nodes_names_out = nodes_genes_intersect + nodes_minus_genes
    expr_df_out = expr_df_in[genes_names_out]
    return(expr_df_out, nodes_names_out, nodes_genes_intersect)



###############################################################################
def rwr_matrix(node_names, network_matrix, restart_matrix, restart_prob, max_iter, tolerance):
    """
    Performs a RWR (Random Walk with Restart) with the given parameters on a
    matrix input. 
    Input: 
        node_names: Name of the nodes in the network
        network_matrix: The probability transition matrix of the network (symmetric)
        restart_matrix: The matrix representing the restart set
        restart_prob: Probability of restart
        max_iter: Maximum number of iterations for convergence
        tolerance: The threshold used with the residual to determine convergence 
    Output:
        num_iter_tmp: Actual number of iterations performed
        residual: The final value of residual
        steady_prob_new: The equlibrium distributions
    """
    no_restart_prob = 1 - restart_prob
    init_prob = 1/len(node_names)
    # Create the vector of probabilities for the nodes
    steady_prob_old = np.empty(np.shape(restart_matrix))
    steady_prob_old.fill(init_prob)
    residual = 100
    num_iter_tmp = 0
    while (residual > tolerance) and (num_iter_tmp < max_iter):
        steady_prob_new = (sparse.csr_matrix.dot(steady_prob_old, network_matrix) \
                            * no_restart_prob + restart_prob \
                            * restart_matrix)
        residual = max(abs(steady_prob_new - steady_prob_old).sum(axis=1))
        print('iteration = ', num_iter_tmp)
        num_iter_tmp += 1
        steady_prob_old = steady_prob_new.copy()
    return(num_iter_tmp, residual, steady_prob_new)

###############################################################################
def rwr_vec(node_names, network_matrix, restart_vec, restart_prob, max_iter, tolerance):
    """
    Performs a RWR (Random Walk with Restart) with the given parameters on a
    vector input. 
    Input: 
        node_names: Name of the nodes in the network
        network_matrix: The probability transition matrix of the network (symmetric)
        restart_vec: The vector representing the restart set
        restart_prob: Probability of restart
        max_iter: Maximum number of iterations for convergence
        tolerance: The threshold used with the residual to determine convergence 
    Output:
        num_iter_tmp: Actual number of iterations performed
        residual: The final value of residual
        steady_prob_new: The equlibrium distribution
    """    # Get the number of nodes
    num_nodes = len(node_names)
    no_restart_prob = 1 - restart_prob
    # Compute the initial probability for the nodes
    init_prob = 1/num_nodes
    # Create the vector of probabilities for the nodes
    steady_prob_old = np.empty(num_nodes)
    steady_prob_old.fill(init_prob)
    # Initialize the loop variables (100 is an arbitrary high value)
    residual = 100
    num_iter_tmp = 0
    while (residual > tolerance) and (num_iter_tmp < max_iter):
        steady_prob_new = sparse.csr_matrix.dot(steady_prob_old, network_matrix)
        steady_prob_new *= no_restart_prob
        steady_prob_new += restart_prob * restart_vec
        # Calculate the residual -- the sum of the absolute
        # differences of the new node probability vector minus the old
        residual = abs(steady_prob_new - steady_prob_old).sum()
        num_iter_tmp += 1
        steady_prob_old = steady_prob_new.copy()
    return (num_iter_tmp, residual, steady_prob_new)

###############################################################################
def import_network(address_net, delimiter):
    """
    Imports the network and generates a dataframe.
    Input:
        address_net: The address of the network
        delimiter: The delimiter used to import the network            
    """
    default_column_headers = ['n_alias_1', 'n_alias_2', 'weight', 'type']
    # Step 1: Read the input
    # the input_file is the network
    with open(address_net, 'r') as fin:
        # Check whether the first line is data or headers
        first_line = fin.readline().strip()
        fin.seek(0)   #go back to the beginning of the file
        fields = first_line.split(sep=delimiter)
        # data
        if is_number(fields[2]):
            net_df = pd.read_csv(fin, sep=delimiter, names=default_column_headers)
        # headers
        else:
            net_df = pd.read_csv(fin, sep=delimiter, header=0)

    # Get the column headers
    node1 = net_df.columns[0]
    node2 = net_df.columns[1]
    weight = net_df.columns[2]
    #t = net_df.columns[3]

    # Get the unique nodes -- the first two columns of the input data,
    # converted to sets to remove duplicates, union'ed, then sorted
    nodes1 = net_df.iloc[:, 0]
    nodes2 = net_df.iloc[:, 1]
    nodes = set(nodes1) | set(nodes2)
    node_names = sorted(nodes)
    num_nodes = len(node_names)
    #print(node_names)

    # Output some info about the input data
    print("Input:")
    print("Lines of data:", len(net_df))
    print("Number of unique nodes:", num_nodes)

    return(node_names, net_df, node1, node2, weight)


###############################################################################
def gen_network_matrix(num_nodes, net_df, node1, node2, weight, node2index):
    """Generates network adjacency matrix and normalizes it"""
    # Transform the first two columns of the DataFrame -- the nodes -- to their indexes
    net_df[node1] = net_df[node1].apply(lambda x: node2index[x])
    net_df[node2] = net_df[node2].apply(lambda x: node2index[x])
    # Create the sparse matrix
    network_matrix = sparse.csr_matrix((net_df[weight].values, (net_df[node1].values, net_df[node2].values)),
                                       shape=(num_nodes, num_nodes), dtype=float)
    # Make the ajdacency matrix symmetric
    network_matrix = (network_matrix + network_matrix.T)
    network_matrix.setdiag(0)
    # Normalize the rows of network_matrix because we are multiplying vector by matrix (from left)
    network_matrix = normalize(network_matrix, norm='l1', axis=1)
    return(net_df, network_matrix)



###############################################################################

def main():
    """The main part of ProGENI"""
    warnings.simplefilter('ignore',SparseEfficiencyWarning)
    ###############################################################################
    args = parse_args()
    
    n_rcg = args.num_RCG
    restart_prob_trans = args.prob_restart_trans
    restart_prob = args.prob_restart_rank
    tolerance = args.tolerance
    max_iter = args.max_iteration
    n_boot = args.num_bootstrap
    percent_boot = args.percent_bootstrap/100
    address_expr = os.path.join(args.directory_expression, args.input_expression)
    address_response = os.path.join(args.directory_response, args.input_response)
    address_net = os.path.join(args.directory_network, args.input_network)
    address_out = os.path.join(args.directory_out, args.output)

    expr_df = pd.read_csv(address_expr, sep=',', header=0, index_col=0).T
    delimiter_net = ','
    (node_names, net_df, node1, node2, weight) = import_network(address_net, delimiter_net)
    
    ###############################################################################
    # Reorder gene column names of expression spreadsheet and gene node names of network
    # such that both name lists start with the intersection of the two list ordered alphabetically
    (expr_df, node_names, nodes_genes_intersect) = spread_match_network(expr_df, node_names)
    
    expr_all = zscore(expr_df.values, axis=0)
    gene_names = list(expr_df.columns.values)
    if len(set(gene_names)) != len(gene_names):
        sys.exit("Duplicate gene names!")
    #    gene2index = {gene_names[i]:i for i in range(len(gene_names))}
    node2index = {node_names[i]:i for i in range(len(node_names))}
    index2node = node_names.copy()
    #    index2gene = gene_names.copy()
    num_nodes = len(node_names)
    
    ###############################################################################
    #Perform network transformation on the gene expression matrix
    (net_df, network_matrix) = \
        gen_network_matrix(num_nodes, net_df, node1, node2, weight, node2index)
    restart_matrix = np.eye(len(node_names), len(nodes_genes_intersect)).T
    print('Obtaining the network')
    (num_iter, residual, gene_similarity_smooth) = \
        rwr_matrix(node_names, network_matrix, restart_matrix, restart_prob_trans, max_iter, tolerance)
    gene_similarity_smooth = gene_similarity_smooth[:, 0:np.size(gene_similarity_smooth, axis=0)]
    gene_similarity_smooth = normalize(gene_similarity_smooth, norm='l1', axis=1)
    expr_all = zscore(expr_all[:, 0:len(nodes_genes_intersect)].dot(gene_similarity_smooth.T), axis=0)
    
    
    ###############################################################################
    ###############################################################################
    response_df = pd.read_csv(address_response, sep=',', header=0, index_col=0).T
    response_all = response_df.values
    names_drug = list(response_df.columns.values)
    n_drug = len(names_drug)
    
    start = time.clock()
    ranked_genes = {}       # Genes ranked such that the most important gene is first
    name_drug_final = []
    aggr_ranked_genes = {}
    # list of interest for drugs:
    #drug_interest = set(names_drug)    #this option keeps all the drugs
    drug_interest = set(['Cisplatin', 'Doxorubicin', '17-AAG'])
    
    
    #######################################################################
    #Obtain global equilibrium distribution of the nodes
    restart_vec = np.ones(len(node_names)) / len(node_names)
    (num_iter, residual, steady_prob_new_global) = \
        rwr_vec(node_names, network_matrix, restart_vec, restart_prob, max_iter, tolerance)
    print(num_iter, residual)
    
    for i_d in range(n_drug):
        if names_drug[i_d] not in drug_interest:
            continue
        y_response = response_all[:, i_d]
        tmp = np.arange(0, len(y_response))
        label_not_nan = list(tmp[~np.isnan(y_response)])
        #we remove samples with missing values
        y_response = y_response[label_not_nan]
        if len(y_response) < 2:
            continue
    
        name_drug_final.append(names_drug[i_d])
        expr_matrix = expr_all[label_not_nan,]
        n_sample = len(label_not_nan)
        n_train = int(round(n_sample * percent_boot))
        ranked_genes[names_drug[i_d]] = []
        np.random.seed(args.seed)
        for k in range(n_boot):
            sample_permute = np.random.permutation(range(n_sample))
            label_train = sorted(sample_permute[np.arange(0, n_train)])
            print('response =', names_drug[i_d], 'repeat =', k)
    
            y_train = y_response[label_train]
            expr_matrix_train = expr_matrix[label_train,]
            #######################################################################
            #Find response-correlated genes (RCG)
            cor_all = [abs(ss.pearsonr(y_train, expr_matrix_train[:, i])[0])
                       for i in range(np.size(expr_matrix_train, 1))]
            argsort_cor = np.argsort(cor_all)
            #######################################################################
            #Use RWR to rank genes w.r.t. the RCG set
            if n_rcg > len(argsort_cor):
                label_feat_seed = argsort_cor.copy()
            else:
                label_feat_seed = argsort_cor[-n_rcg:]  #Index of most correlated features
            restart_vec = np.zeros(len(node_names))
            for i in range(len(label_feat_seed)):
                if gene_names[label_feat_seed[i]] in index2node:
                    restart_vec[node2index[gene_names[label_feat_seed[i]]]] = \
                        cor_all[label_feat_seed[i]]
            restart_vec /= restart_vec.sum()
    
            (num_iter, residual, steady_prob_new_drug) = \
                rwr_vec(node_names, network_matrix, restart_vec,
                        restart_prob, max_iter, tolerance)
            print('RWR iterations =', num_iter, 'Residual = ', residual)
            steady_prob_new = steady_prob_new_drug - steady_prob_new_global
            ranked_nodes = np.argsort(steady_prob_new)[::-1]     #indices sorted high prob to low
            #Ranked gene names with the same order as above line except that
            #it removes genes not in the gene expression dataset
            ranked_names = [index2node[i] for i in ranked_nodes if index2node[i] in gene_names]
    
            ranked_genes[names_drug[i_d]].append(ranked_names)
        
        aggr_ranked_genes[names_drug[i_d]] = \
            rank_aggregate_borda(ranked_genes[names_drug[i_d]], 'geometric_mean')
        aggr_ranked_genes_df = pd.DataFrame(data=aggr_ranked_genes, index=None, columns=drug_interest)
        aggr_ranked_genes_df.to_csv(address_out, sep=',')

    end_alg = time.clock()
    print(end_alg-start)



if __name__ == "__main__":
    main()
