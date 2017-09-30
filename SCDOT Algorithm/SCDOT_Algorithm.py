import pandas as pd
import matplotlib.pyplot as plt
from graphtools import *
from collections import Counter
# from visualization import *
# from chameleon import *
# import scipy.io as scio
# import seaborn as sns
import sys
sys.path.append('../')
from NMI import *
if __name__ == "__main__":
    # get a set of data points
    # data = scio.loadmat('./datasets/can383.mat',header=None, names=['x', 'y'])
    df = pd.read_csv('./datasets/Wholescale.csv', delimiter=' ', header=None)
    # create knn graph
    k=30;cluster_number=2
    rho, dismatrix,k_neighbors = Rho_Calculate(df, k)
    Delta,nneigh=Delta_Calculate(rho, dismatrix)
    graph=Creat_graph(Delta,nneigh,df)


    L_subgraph,L_outlier=part_graph_1(graph)
    print(len(L_subgraph))
    print(len(L_outlier))
    Sub_clu_num=len(L_subgraph)
    df['cluster'] = 0;cn=1;flag=0;ind=1
    num_subgraph=range(Sub_clu_num)

    con=connectivity_cal(L_subgraph,k_neighbors,dismatrix)
    while True:

        con,df,cn,flag,L_subgraph,num_subgraph=merge_graph(con,df,L_subgraph,cn,num_subgraph)
        if flag==1 or len(Counter(num_subgraph))==cluster_number :
            break

    print(len(Counter(num_subgraph)))
    num=Counter(num_subgraph).items()
    target=[i[0] for i in num if i[1]==1]
    for i in range(len(target)):
        g=target[i]
        li=L_subgraph[g].nodes()
        df.ix[li, 'cluster'] = cn
        cn=cn+1
    L_outlier=[i.nodes() for i in L_outlier]
    L_outlier=[i[0] for i in L_outlier]
    df.ix[L_outlier, 'cluster']=cn
    temp=df['cluster']
    sort_class=Counter(temp).items()
    for i in range(len(sort_class)):
        index=[a for a in range(len(df)) if df.ix[a,'cluster'] == sort_class[i][0]]
        df.ix[index, 'cluster']=ind
        ind=ind+1
    # df.plot(kind='scatter', c=df['cluster'], cmap='gist_rainbow', x='x', y='y')
    label=pd.read_csv('./datasets/Wholescale.csv',delimiter=' ',header=None)
    print(NMI(np.array(label[0]),df['cluster'].values.T))
    plt.show()
