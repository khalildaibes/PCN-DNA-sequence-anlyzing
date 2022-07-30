import pandas as pd
from tqdm import tqdm

tqdm.pandas()
import operator


def secondG(dff, name):
    print("second clustering starts...")
    maxx = dff.ClusterID.max()

    def prkmer(list1, j):
        AA = {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0,
              'F': 0, 'P': 0, 'S': 0, 's': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}
        new = ''
        for i in range(20):
            for k in list1:
                if k[i] != 'x':
                    AA[k[i]] += 1
            new += max(AA.items(), key=operator.itemgetter(1))[0]
            AA = dict.fromkeys(AA, 0)
        new += ':'
        new += str(j)
        return new

    lol = [dff[dff['ClusterID'] == i].kmer.tolist() for i in tqdm(range(maxx + 1))]
    newl = [[prkmer(lol[i], i), i] for i in tqdm(range(len(lol)))]

    newdf = pd.DataFrame(newl, columns=['kmer', 'id'])
    klis = newdf.kmer.tolist()
    klis = list(set(klis))
    st1 = {}
    for l in tqdm(klis):
        st1[l] = 0

    def hamming(a, b):
        cnt = 0
        for i in range(20):
            if a[i] == b[i] and a[i] != 'x' and b[i] != 'x':
                cnt += 1
                if cnt > 10:
                    return 1
        return 0

    cl = -1
    data = []
    for i in tqdm(range(len(klis))):
        if st1[klis[i]] != 0:
            continue
        cl += 1
        data.append([klis[i], cl])
        for j in range(i + 1, len(klis)):
            if st1[klis[j]] == 0:
                if hamming(klis[i], klis[j]) == 1:
                    st1[klis[j]] += 1
                    data.append([klis[j], cl])

    dat = pd.DataFrame(data, columns=['kmer', 'clusterid'])
    print(dat.head())

    dat.to_csv(name, index=False)



