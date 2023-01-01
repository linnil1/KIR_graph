"""linnil1 modified from 02_determine_copy_number.py"""
from scipy.signal import argrelextrema
from sklearn.neighbors import KernelDensity
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def calcPloidyThreshold(cov, name):
    colors = ['r','g','b','c','m','y']
    thres_dic = {}
    for i in range(len(cov)):
        gene = cov.index[i]
        a = np.array(cov.iloc[i,:]).reshape(-1,1)
        kde = KernelDensity(kernel='gaussian', bandwidth=0.075).fit(a)
        s = np.linspace(min(a)-0.05,max(a)+0.05)
        e = kde.score_samples(s.reshape(-1,1))
        mi, ma = argrelextrema(e,np.less)[0], argrelextrema(e,np.greater)[0]
        n_thress = len(mi)
        n_bin = len(ma)
        thres_dic[gene] = s[mi]
        if False and n_thress > 0:
            plt.figure()
            plt.plot(s[:mi[0]+1],e[:mi[0]+1], colors[0])
            for j in range(1,n_bin-1):
                plt.plot(s[mi[j-1]:mi[j]+1],e[mi[j-1]:mi[j]+1], colors[j])
            plt.plot(s[mi[n_bin-2]:],e[mi[n_bin-2]:],colors[(n_bin-1) % len(colors)])
            plt.savefig(f"{name}.{gene}.png")
            print(f"Save {name}.{gene}.png")
        else:
            print(gene + ' had zero threshold')
    return thres_dic


def calcPloidy(cov, thres_dic):
    # 02_determine_copy_number.py
    genelist = ['KIR3DS1', 'KIR3DL1', 'KIR2DS4', 'KIR2DS3;KIR2DS5', 'KIR2DS2',
                'KIR2DS1', 'KIR2DP1', 'KIR2DL5A;KIR2DL5B', 'KIR2DL3', 'KIR2DL2',
                'KIR2DL1', 'KIR3DL3', 'KIR3DL2', 'KIR2DL4']
    tmp = np.zeros((len(genelist),len(cov.columns)),dtype = "int")
    copy = pd.DataFrame(tmp)
    copy.columns = cov.columns
    copy.index = genelist
    for gene in genelist:
        cut_thres = np.hstack(([0], np.ravel(thres_dic[gene]), [4]))
        ratio = cov.loc[gene, :]
        # the original method may fail
        # if gene != "KIR2DL4":
        #     copy.loc[gene, :] = np.array(pd.cut(ratio, cut_thres, labels=False))
        # else:
        #     copy.loc[gene, :] = np.array(pd.cut(ratio, cut_thres, labels=np.array([1,2,3,4])))
        copy.loc[gene, :] = np.array(pd.cut(ratio, cut_thres, labels=False))
    return copy


def getPloidy(cov: pd.DataFrame, name: str) -> pd.DataFrame:
    thres_dic = calcPloidyThreshold(cov, name)  # type: ignore
    return calcPloidy(cov, thres_dic)  # type: ignore
