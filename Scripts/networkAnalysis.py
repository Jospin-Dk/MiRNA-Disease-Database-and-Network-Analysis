##################################################################################################
#PACKAGES
##################################################################################################
import pandas as pd
import networkx as nx

##################################################################################################
#DATA
##################################################################################################
print("####################")
print("LOADING DATABASE")
print("####################")

#reading the data
miRNA_locus_gene_disease=pd.read_csv("./miRNA_gene_disease.csv")
#removing 1st column (Unnamed: 0)
miRNA_locus_gene_disease.drop('Unnamed: 0', inplace=True, axis=1)

print("####################")
print("DATABASE LOADED")
print("####################")

##################################################################################################
#NETWORK BUILDING
##################################################################################################
#instantiating a graph G
G= nx.Graph()

#extracting unique miRNA, genes , diseases
unique_mirna= pd.unique(miRNA_locus_gene_disease['miRNA'])
unique_gene= pd.unique(miRNA_locus_gene_disease['geneSymbol'])
unique_disease= pd.unique(miRNA_locus_gene_disease['diseaseName'])

#adding all unique miRNA, genes, diseases as nodes in the graph
G.add_nodes_from(unique_mirna)
G.add_nodes_from(unique_gene)
G.add_nodes_from(unique_disease)

#assigning edges
for x in range(0,miRNA_locus_gene_disease.shape[0]):
    G.add_edges_from([(miRNA_locus_gene_disease.iloc[x,0],miRNA_locus_gene_disease.iloc[x,4])]) #adding an edge between each miRNA(0) and gene(4)
    G.add_edges_from([(miRNA_locus_gene_disease.iloc[x,4],miRNA_locus_gene_disease.iloc[x,5])]) #adding an edge between each gene(4) and disease(5)

##################################################################################################
#MAIN
##################################################################################################
#answer for question 1
print("Question 1:")
print(nx.info(G))
print("----------")

#answer for question 2
print("Question 2:")
bridges=list(nx.bridges(G))
print("The Graph has", len(bridges),"gene connections that are bridges")
print("----------")

#answer for question 3
print("Question 3: ")
most_connected=max(dict(G.degree()).items(), key=lambda x : x[1])
print("The most connected node (hub) is:", most_connected[0])
print("The hub has:",most_connected[1],"connections")
print("----------")

#answer for question 4
print("Question 4:")
betw_cent= nx.betweenness_centrality(G)
high_betw_cent=max(betw_cent, key=betw_cent.get)
print("The most important bottleneck node (highest betweenness centrality) is: ", high_betw_cent)
print("----------")

#answer for question 5
print("Question 5:")
close_cent= nx.closeness_centrality(G)
high_close_cent= max(close_cent, key=close_cent.get)
print("The most influential node (highest closeness centrality) is: ", high_close_cent)
print("----------")


