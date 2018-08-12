#usage python dcor.py data 1 4 4 7

import sys
from random import sample
import numpy as np

data_file=open(sys.argv[1],"r").readlines()
gene_1_index=range((int(sys.argv[2])),(int(sys.argv[3])))
gene_2_index=range((int(sys.argv[4])),(int(sys.argv[5])))
case_set=[]
control_set=[]
for l in data_file:
    l=[int(x) for x in l.split(" ")]
    if l[0]==1:
        case_set+=[l]
    else:
        control_set+=[l]

def dist(x, y, idx):
    return sum([(x[i]-y[i])**2 for i in idx])
        
def dcor(dataset, idx1, idx2):
    k=len(dataset)
    mat1=[]
    mat2=[]
    for i in range(k):
        mat1+=[[]]
        mat2+=[[]]
        for j in range(k):
            mat1[i]+=[dist(dataset[i], dataset[j], idx1)]
            mat2[i]+=[dist(dataset[i], dataset[j], idx2)]
    sumrow1=[0]*k
    sumrow2=[0]*k
    sumcol1=[0]*k
    sumcol2=[0]*k
    sumall1=0
    sumall2=0
    for i in range(k):
        for j in range(k):
            sumrow1[i]+=mat1[i][j]
            sumrow2[i]+=mat2[i][j]
            sumcol1[j]+=mat1[i][j]
            sumcol2[j]+=mat2[i][j]
            sumall1+=mat1[i][j]
            sumall2+=mat2[i][j]
    for i in range(k):
        for j in range(k):
            mat1[i][j]-=(sumrow1[i]*1.0/k+sumcol1[j]*1.0/k-sumall1*1.0/k/k)
            mat2[i][j]-=(sumrow2[i]*1.0/k+sumcol2[j]*1.0/k-sumall2*1.0/k/k)
    cov12=0
    cov1=0
    cov2=0
    for i in range(k):
        for j in range(k):
            cov1+=mat1[i][j]**2
            cov2+=mat2[i][j]**2
            cov12+=mat1[i][j]*mat2[i][j]
    #print cov1, cov2, cov12
    return cov12/((cov1*cov2)**0.5)

n1=len(case_set)
n2=len(control_set)
k=50
m=10
nz=10
nmc=5

def U_dcor(m, k, data_set, idx1, idx2):
    U=0
    for i in range(m):
        subsample=sample(data_set, k)
        U+=dcor(subsample, idx1, idx2)
    return U/m

def estimate_zetas(nz, nmc, k, data_set, idx1, idx2):
    res1=np.array([0.]*nz)
    res2=np.array([0.]*nz)
    for i in range(nz):
        subsample=sample(data_set, 1)
        for j in range(nmc):
            subsample+=sample(data_set, k-1)
            res=dcor(subsample, idx1, idx2)
            #print i, j, res
            if j==0:
                res1[i]=res
            res2[i]=res2[i]+res/nmc
    return (np.cov(res1), np.cov(res2))

zetas_case=estimate_zetas(nz, nmc, k, case_set, gene_1_index, gene_2_index)
zetas_control=estimate_zetas(nz, nmc, k, control_set, gene_1_index, gene_2_index)

var_case=(k**2*m*zetas_case[1]/n1+zetas_case[0])/m
var_control=(k**2*m*zetas_control[1]/n2+zetas_control[0])/m

s1=U_dcor(m, k, case_set, gene_1_index, gene_2_index)
s2=U_dcor(m, k, control_set, gene_1_index, gene_2_index)

print s1, s2, var_case, var_control, "t=", (s1-s2)/(var_case+var_control)**0.5
