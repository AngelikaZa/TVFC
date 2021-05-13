## Library of useful functions for connectomics
# Connectome Library 

##### Required modules
# Import all necessary libraries
import pandas as pd
import numpy as np
import scipy.stats as stats
import brainconn as con
import bct as bct
import networkx as nx


def calculate_edges(list_of_paths):
	'''
	Aim: A function to calculate the average number of edges in a population of graphs
	Inputs: 
		list_of paths: a list of paths for the individual participant connectomes
	Outputs:
		float: average number of connections (edges)
	'''
	edges = []
	for i in range(len(List)):
		path = list_of_paths[i]
		data = pd.read_csv(path, header=None, sep=" ")
		x = np.count_nonzero(data.values)
		edges.append(x)
		num = (np.mean(np.array(edges))) / 2
	return(num)


def module_detection(path_averageConnectome, g, iterations):
	''' 
	Aim: A function to identify network modules using Louvain community algorithm and consensus approach after 1000 iterations
	Inputs: 
		path_averageConnectome: the path to the average connectome for all participants
		g: desired gamma (float)
		iterations: number of desired iterations for the Louvain algorithm
	Outputs: 
		1. connectivity matrix
		2. vector allocating each node to a module
	Details: G values: small modules (<1) or large modules (>1), default = 1
 	''' 
	vectors = []
	qstats = []
	# Run Louvain algorithm
	for i in range(iterations):
		df = pd.read_csv(path_averageConnectome, header=None) 
		vector, qstat = con.modularity.community_louvain(df.values, gamma=g)
		vectors.append(vector)
		qstats.append(qstat)
	# Construct agreement matrix
	qstats = np.array(qstats)
	vectors = np.array(vectors)
	vectors = np.moveaxis(vectors, 0, -1)
	agr_matrix = con.clustering.agreement(vectors, buffsz=150)
	agr_matrix = con.utils.normalize(agr_matrix)
	# Consensus approach
	cons_matrix = con.clustering.consensus_und(agr_matrix, 0.5, 100)
	return(cons_matrix, vectors)


def module_connections(graph, module1, module2):
	'''
	Aim: A function to calculate the sum strength of connection between modules
	Inputs
		graph: as a pandas dataframe
		2 modules (each a list of the respective nodes allocated to each module)  	
	Outputs: 
		sum strength (integer) 
	'''
	graph = pd.DataFrame(graph, index=module1)
	graph = graph[module2]
	sumStrength = graph.values.sum()
	return sumStrength


def module_values(graph, module1, module2): 
	'''
	Aim: A function to extract a list of total values from a subselection of a matrix
	Inputs: 
		graph (as a pandas dataframe) and 2 modules
	Outputs: l
		list of connection values from the subselection of the graph
	'''
	graph = pd.DataFrame(graph, index=module1)
	graph = graph[module2]
	values = graph.values
	values = values.flatten()
	return values

def module_length(graph, module1, module2):
	'''
	Aim: A function to extract average characteristic path length from a subsection of a matrix
	Inputs: a graph and 2 modules
	Outputs: mean connection length for that module (float)
	'''
	graph = pd.DataFrame(graph, index=module1)
	graph = graph[module2]
	meanLength = graph.values.sum() / (len(module1, module2))
	return meanLength



def tahn_transform(graph, average_connectome, std_connectome):
	'''
	Aim: function to z and tahn transform the values of a connection to atrophy values.
	Inputs: 
		graph: connectivity strength connectome (per individual)
		average connectome: a control average connectome
		std_connectome: control std connectome. 
		All connectomes need to be in pandas dataframe format.
	Outputs: 
		connectivity matrix with tahn transformed values inplace of connectivity strength
	'''
	data = np.genfromtxt(graph)
	data_z = (average_connectome - data) / std_connectome
	data_transform = np.tahn(data_z)
	return data_transform


def annotation(Path):#
	'''
	Aim: Loading and extracting data from a matlab file containing gene expression data from Allen atlas in Glasser space (extracted using matlab code from Arnatkevic̆i et al. Neuroimage, 2019)
	Inputs: the path to the matlab file
	Outputs: gene list, gene symbol list
	''' 
	annot = sio.loadmat(Path)
	names = annot["probeInformation"] #find probe information on dictionary 
	df = pd.DataFrame(data=names.flatten()) #load as dataframe
	geneList = df["EntrezID"] #find gene symbol list
	geneList= geneList[0].tolist() # transform dataframe column to list
	geneSymbol= df["GeneSymbol"]
	geneSymbol = geneSymbol[0].tolist()
	return geneList, geneSymbol


def spinPermutation(permutation, measure1, measure2, n_permutations = 1000):
	'''
	Aim: Calculate p-spin between 2 imaging measures measure 1 and measure 2 compared to spatial nulls
	Inputs: 
		permutation: np.array of permuted indices - output of spin permutation done in R
		measure1: np.array of parcellated imaging measure
		measure2: np.array of parcellated imaging measure
		n_permutations: the number of permutations to run, default = 1000
	Outputs: 
		tuple (r, p-spin): 
			r: the spearman correlation coefficient between measure 1 and measure 2 
			p-spin: the p value 
	'''
	n_roi = len(measure1)
	r_spin = np.empty(n_permutations)
	r_obs, p_obs = stats.spearmanr(measure1, measure2)

	for i in range(n_permutations):
		rotated_index = permutation[i]
		perm = np.empty(n_roi)
		for j in range(n_roi):
			perm[j] = measure1[(rotated_index[j] - 1)] # -1 as starting from 0
			r_spin[i] = stats.spearmanr(perm, measure2)[0]
	pv_spin = np.mean(np.abs(r_spin) >= np.abs(r_obs))
	return(r_obs, pv_spin)

def calculateDensity(matrix):
	'''
	Aim: calculate the density of a connectome
	Inputs:
		connectivity matrix as np array
	Outputs: 
		density of the matrix (float)
	'''
	n = len(matrix)
	dens_net = sum(sum(matrix)) / (matrix.max(0) * n * (n-1)).max()
	return(dens_net)


def calculateEntropy(matrix):
	'''
	Aim: calculate the entrophy of a connectome
	Inputs: 
		conectivity matrix as pd.dataframe
	Outputs: 
		the mean entropy (float)
	'''
  	## First bin each node's connectivity pattern
	n = len(matrix)
	binned = pd.DataFrame(np.zeros((n,n)))
	entropy = []
	# calculate shannon entropy for each node
	for i in range(len(matrix)):  
		binned[i] = pd.cut(matrix[i], bins=10, labels=False)
	# normalise compared to entropy of normal distribution
	ent = scipy.stats.entropy(binned[i], base=10) / np.log10(n)
	entropy.append(ent)
	# mean entropy for the node
	mean_entropy = np.mean(entropy)
	return(mean_entropy)