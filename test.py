# coding: utf-8

# ipython notebook requires this
# %matplotlib inline

# python console requires this
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import venn
q_val = ["01","001","1","005","05"]
for q in q_val:
    reg = [1000,2000,3000,4000,5000]
    dfs_genes = []
    for r in reg:
        fname = "intersect/node_list_"+str(r)+"_"+q+".tsv"
        df = pd.read_csv(fname,sep="\t")
        df_genes = set(df.loc[df["Type"]=="gene","Label"])
        dfs_genes.append(df_genes)
    labels = venn.get_labels(dfs_genes, fill=['number', 'logic'])
    fig, ax = venn.venn5(labels, names=reg)
    fig.savefig('venn'+q+'.png', bbox_inches='tight')
