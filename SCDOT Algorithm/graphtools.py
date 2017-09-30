import numpy as np
import networkx as nx
from tqdm import tqdm

import metis

def euclidean_distance(a, b):
    return np.linalg.norm(np.array(a) - np.array(b))

def knn_graph(df, k, verbose=False):
    points = [p[1:] for p in df.itertuples()]
    g = nx.Graph()
    if verbose: print "Building kNN graph (k = %d)" % (k)
    iterpoints = tqdm(enumerate(points), total=len(points)) if verbose else enumerate(points)
    for i, p in iterpoints:
        distances = map(lambda x: euclidean_distance(p, x), points)
        closests = np.argsort(distances)[1:k+1] # second trough kth closest
        for c in closests:
            g.add_edge(i, c, weight=distances[c])
        g.node[i]['pos'] = p
    return g
def Rho_Calculate(df,k):
    r,c=df.shape
    points = [p[1:] for p in df.itertuples()]
    iterpoints = tqdm(enumerate(points), total=len(points))
    dismatrix = np.zeros([r, r])
    rho = np.zeros(r)
    k_neighbors = [[0] * k for i in range(r)]
    for i, p in iterpoints:
        sum=0.0

        distances = map(lambda x: euclidean_distance(p, x), points)
        ord_pre = [(a, b) for a, b in enumerate(distances)]
        ord = sorted(ord_pre, key=lambda x: x[1],reverse=False)
        k_neighbors[i][:]=[j[0] for j in ord[1:k+1]]
        closests = np.argsort(distances)[1:k+1] # second trough kth closest

        for c in closests:
            sum=sum+distances[c]
        rho[i]=1/sum
        dismatrix[i,:]=distances
    return rho,dismatrix,k_neighbors
def Delta_Calculate(rho,dismatrix):
    maxd=np.amax(dismatrix)
    r,c=dismatrix.shape
    r=int(r)
    Delta=np.zeros(r);nneigh=np.zeros(r)
    ord_pre=[(a,b) for a,b in enumerate(rho)]
    ord=sorted(ord_pre,key=lambda x:x[1],reverse=True)
    if rho[ord[0][0]]==rho[ord[1][0]]:
        rho[ord[0][0]]=rho[ord[0][0]]+0.0001
    Delta[ord[0][0]]=0
    nneigh[ord[0][0]]=-1
    for ii in range(1,r):
        Delta[ord[ii][0]]=maxd
        for jj in range(0,ii):
            if dismatrix[ord[ii][0]][ord[jj][0]]<Delta[ord[ii][0]]:
                Delta[ord[ii][0]]=dismatrix[ord[ii][0]][ord[jj][0]]
                nneigh[ord[ii][0]]=ord[jj][0]
    Delta[ord[0][0]] = np.amin(dismatrix[ord[0][0]][:]);
    nneigh[ord[0][0]] = ord[0][0]
    nneigh=[int(i) for i in nneigh]
    return   Delta,nneigh

def Creat_graph(Delta,nneigh,df):
    m=len(nneigh)
    g=nx.Graph()
    for i in range(m):
        g.add_edge(i,nneigh[i],weight=Delta[i])
        g.node[i]['pos']=tuple(df.ix[i])
    return g
def box_method(graph):
    n=graph.number_of_nodes()
    list=[]
    edge_inf=graph.edges()
    for i in range(n):
        list.append(graph.get_edge_data(edge_inf[i][0],edge_inf[i][1])['weight'])
    list.sort()
    Q1=list[int(n/4)]
    Q3=list[int(n*3/4)]
    threshold_1=Q3+1.5*(Q3-Q1)
    threshold_2=Q3+3*(Q3-Q1)
    return threshold_1,threshold_2
def part_graph_1(graph):
    th1,th2=box_method(graph)
    n = graph.number_of_nodes()
    edge_inf = graph.edges()
    for i in range(n):
        if graph.get_edge_data(edge_inf[i][0],edge_inf[i][1])['weight']>th2:
            graph.remove_edge(edge_inf[i][0],edge_inf[i][1])
    g1 = nx.connected_component_subgraphs(graph)
    L=[]
    L_outlier=[]
    while True:
        try:
            sub_graph = next(g1)
            # L=part_graph_1(sub_graph,L)
            if sub_graph.number_of_nodes()>1:
                L.append(sub_graph)
            if sub_graph.number_of_nodes()==1:
                L_outlier.append(sub_graph)
        except StopIteration:
            break
    return L,L_outlier
def connectivity_cal(L,k_neighbors,dismatrix):
    n=len(L)
    con=[[0]*n for i in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j:
                con[i][j]=0
                continue
            con[i][j]=subgraph_connec_calculate(L[i],L[j],k_neighbors,dismatrix)
        con[i][i]=0
    return con
def subgraph_connec_calculate(graph1,graph2,k_neighbors,dismatrix):
    l1=graph1.nodes()
    l2=graph2.nodes()
    b1=0;b2=0;sum=0
    for i in l1:
        for j in l2:
            if i in k_neighbors[j]:
                b1=1
            if j in k_neighbors[i]:
                b2=1
            sum=sum+(b1+b2)/dismatrix[i][j]
    dis=sum/len(l1)/len(l2)
    return dis


    # return
def part_graph_2(graph):
    th1,th2=box_method(graph)
    n = graph.number_of_nodes()
    edge_inf = graph.edges()
    for i in range(n):
        if graph.get_edge_data(edge_inf[i][0],edge_inf[i][1])['weight']>th1:
            graph.remove_edge(edge_inf[i][0],edge_inf[i][1])
    g1 = nx.connected_component_subgraphs(graph)
    L=[]
    while True:
        try:
            sub_graph = next(g1)
            if sub_graph.number_of_nodes()>1:
                L.append(sub_graph)
        except StopIteration:
            break
    return L
# def part_graph_3(graph):
def part_graph(graph, k, df=None):
    edgecuts, parts = metis.part_graph(graph, k)
    for i, p in enumerate(graph.nodes()):
        graph.node[p]['cluster'] = parts[i]
    if df is not None:
        df['cluster'] = nx.get_node_attributes(graph, 'cluster').values()
    return graph

def get_cluster(graph, clusters):
    nodes = [n for n in graph.node if graph.node[n]['cluster'] in clusters]
    return nodes

def connecting_edges(partitions, graph):
    cut_set = []
    for a in partitions[0]:
        for b in partitions[1]:
            if a in graph:
                if b in graph[a]:
                    cut_set.append((a, b))
    return cut_set

def min_cut_bisector(graph):
    graph = graph.copy()
    graph = part_graph(graph, 2)
    partitions = get_cluster(graph, [0]), get_cluster(graph, [1])
    return connecting_edges(partitions, graph)
def merge_graph(con,df,L,cn,num_subgraph):
    temp=0;flag=0
    for i in range(len(con)):
        for j in range(len(con)):
            if con[i][j]>temp:
                temp=con[i][j]
                a,b=i,j
    if temp>0:
        graph1,graph2=L[a],L[b]
        con[a][b]=0;con[b][a]=0
        l1=graph1.nodes();l2=graph2.nodes()
        l=l1+l2
        if df.ix[l[0],'cluster']!=0 or df.ix[l[-1],'cluster']!=0:
            if  df.ix[l[0],'cluster']<df.ix[l[-1],'cluster']:
                if df.ix[l[0],'cluster']==0:
                    df.ix[l, 'cluster'] = df.ix[l[-1],'cluster']
                else:
                    df.ix[l, 'cluster'] = df.ix[l[0], 'cluster']
            else:
                if df.ix[l[-1],'cluster']==0:
                    df.ix[l, 'cluster'] = df.ix[l[0],'cluster']
                else:
                    df.ix[l, 'cluster'] = df.ix[l[-1], 'cluster']
        else:
            df.ix[l, 'cluster']=cn
            cn=cn+1
        graph1.add_path(graph2.nodes())
        L[a] = graph1
        L[b] = graph1
        index=[i for i, v in enumerate(num_subgraph) if v ==num_subgraph[b]]
        for j in index:
            num_subgraph[j] = num_subgraph[a]

            # L.pop(b)
    else:
        flag=1


    return con,df,cn,flag,L,num_subgraph
def get_weights(graph, edges):
    return [graph[edge[0]][edge[1]]['weight'] for edge in edges]

def bisection_weights(graph, cluster):
    cluster = graph.subgraph(cluster)
    edges = min_cut_bisector(cluster)
    weights = get_weights(cluster, edges)
    return weights