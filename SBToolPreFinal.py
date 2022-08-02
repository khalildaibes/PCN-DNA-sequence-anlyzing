import pandas as pd
from tqdm import tqdm

tqdm.pandas()


def indx(row):
    row = row[21:]
    return row


def wrap_indx(x):
    return indx(x['kmer'])


def linkdata(df, ref1, ref2, name):
    print("start linking data...")
    selected_columns = df[["kmer"]]
    df1 = selected_columns.copy()
    df['kmerid'] = df1.progress_apply(wrap_indx, axis=1)

    maxx = df.clusterid.max()

    bigc = []
    for i in tqdm(range(maxx + 1)):
        newl = []
        flat_list = []
        d = df[df['clusterid'] == i]
        blist = d.kmerid.tolist()
        for i in blist:
            a = ref2.loc[ref2['ClusterID'] == int(i)]
            a = a.id.tolist()
            newl.append(a)
        flat_list = [item for sublist in newl for item in sublist]
        bigc.append(flat_list)

    big = pd.DataFrame()
    for i in tqdm(range(len(bigc))):
        aa = ref1.iloc[bigc[i]]
        l = aa.kmer.tolist()
        aa.insert(7, 'cluster', str(i))
        aa = aa.reset_index(drop=True)
        big = pd.concat([big, aa.copy()])

    x = big.cluster.value_counts()
    x = x.sort_index()
    p = pd.DataFrame(x)
    p = p.reset_index()
    p = p.rename(columns={'index': 'clnum', 'cluster': 'count'}, inplace=False)
    lastlist = p.values.tolist()
    num = {}
    for l in tqdm(lastlist):
        num[l[0]] = l[1]

    def nums(x):
        return num[x]

    def pos(x):
        x = x[1:4]
        if x[2] == ',':
            return int(x[:2])
        else:
            return int(x)

    big['cnt'] = big.apply(lambda row: nums(row['cluster']), axis=1)
    big['pos'] = big.apply(lambda row: pos(row['position']), axis=1)
    print("number of final kmers:" + str(len(big)))
    print(big.head())
    big.to_csv(name, index=False)




