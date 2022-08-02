import pandas as pd
from tqdm import tqdm

tqdm.pandas()


def firstG(kmers, name):
    print("first clustering starts...")
    selected_columns = kmers[["kmer"]]
    df = selected_columns.copy()

    ## functions
    def fitri(row):
        row = row[:3]
        return row

    def latri(row):
        row = row[17:]
        return row

    def betw(row):
        row = row[3:17]
        return row

    def wrap_f(x):
        return fitri(x['kmer'])

    def wrap_l(x):
        return latri(x['kmer'])

    def wrap_b(x):
        return betw(x['kmer'])

    def diff_letters(a, b):
        cnt = 0
        for i in range(len(a)):
            if a[i] == b[i] and a[i] != 'x' and b[i] != 'x':
                cnt += 1
                if cnt > 4:
                    return 1
        return 0

    df['first'] = df.progress_apply(wrap_f, axis=1)
    df['last'] = df.progress_apply(wrap_l, axis=1)
    df['k'] = df.index
    df['kmers'] = df.progress_apply(wrap_b, axis=1)
    df['id'] = kmers['id'].copy()

    lastlist = df.values.tolist()
    ids = {}
    for l in tqdm(lastlist):
        ids[l[3]] = l[5]

    st = {}
    for l in tqdm(lastlist):
        st[l[3]] = 0

    km = {}
    for l in tqdm(lastlist):
        km[l[3]] = l[0]

    k = {}
    for l in tqdm(lastlist):
        k[l[3]] = l[4]

    last = {}
    for l in tqdm(lastlist):
        last[l[3]] = l[2]

    first = {}
    for l in tqdm(lastlist):
        first[l[3]] = l[1]

    firstd = {}
    for i in tqdm(lastlist):
        firstd[i[1]] = []

    for j in tqdm(lastlist):
        firstd[j[1]].append(j[3])

    lastd = {}
    for i in tqdm(lastlist):
        lastd[i[2]] = []

    for j in tqdm(lastlist):
        lastd[j[2]].append(j[3])

    cl = -1
    data = []
    for l in tqdm(firstd.values()):
        for i in range(len(l)):
            if st[l[i]] != 0:
                continue
            cl += 1
            data.append([km[l[i]], cl, ids[l[i]]])
            for j in range(i + 1, len(l)):
                if last[l[i]] == last[l[j]]:
                    if diff_letters(k[l[i]], k[l[j]]) == 1:
                        st[l[j]] += 1
                        data.append([km[l[j]], cl, ids[l[j]]])
                    else:
                        st[l[j]] += 1
                        data.append([km[l[j]], -1, ids[l[j]]])

    dff = pd.DataFrame(data, columns=['kmer', 'ClusterID', 'id'])

    dff = dff[dff['ClusterID'] != -1]
    dff = dff.reset_index(drop=True)

    print("number of kmers after first clusterng: " + str(len(dff)))
    print(dff.head())
    dff.to_csv(name, index=False)




