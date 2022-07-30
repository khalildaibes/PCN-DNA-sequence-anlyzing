import pandas as pd
from tqdm import tqdm

tqdm.pandas()


def diff_letters(a, b):
    cnt = 0
    for i in range(len(a)):
        if a[i] != b[i] and a[i] != 'x' and b[i] != 'x':
            cnt += 1
            if cnt > 1:
                break
    return cnt


def fiten(row):
    row = row[:10]
    return row


def laten(row):
    row = row[10:]
    return row


def wrap_f(x):
    return fiten(x['kmer'])


def wrap_l(x):
    return laten(x['kmer'])


def var(kmers0, kmers1, name):
    print("Third step (1-mismatch combined k-mers):\n")
    kmers = kmers1.copy()
    kmers1['tmpid'] = kmers1.index
    kmers1['var'] = 0

    print("number of k-mers before combination: " + str(len(kmers0)))

    kmers1['id'] = kmers0['id']

    selected_columns = kmers1[["tmpid", "var"]]
    new_df = selected_columns.copy()

    all_list = new_df.values.tolist()
    All = {}
    st = {}
    for l in tqdm(all_list):
        All[l[0]] = l[1]
        st[l[0]] = l[1]

    kmers['first'] = kmers.progress_apply(lambda r: fiten(r['kmer']), axis=1)
    kmers['last'] = kmers.progress_apply(lambda r: laten(r['kmer']), axis=1)
    kmers['kmer'] = kmers.index

    lastlist = kmers.values.tolist()
    removed = {}
    for i in tqdm(lastlist):
        removed[i[0]] = []
    last = {}
    for l in tqdm(lastlist):
        last[l[0]] = l[2]
    first = {}
    for l in tqdm(lastlist):
        first[l[0]] = l[1]
    firstd = {}
    for i in tqdm(lastlist):
        firstd[i[1]] = []
    for j in tqdm(lastlist):
        firstd[j[1]].append(j[0])
    lastd = {}
    for i in tqdm(lastlist):
        lastd[i[2]] = []
    for j in tqdm(lastlist):
        lastd[j[2]].append(j[0])

    def mapf():
        r = []
        for l in tqdm(firstd.values()):
            if len(l) > 1:
                for i in range(len(l)):
                    if st[l[i]] != 0:
                        continue
                    for j in range(i + 1, len(l)):
                        if st[l[j]] == 0:
                            if diff_letters(last[l[i]], last[l[j]]) <= 1:
                                All[l[i]] += 1
                                st[l[j]] += 1
                                removed[l[i]].append(l[j])
                                r.append(l[j])

        # second comparision
        for e in tqdm(lastd.values()):
            if len(e) > 1:
                for x in range(len(e)):
                    if st[e[x]] != 0:
                        continue
                    for y in range(x + 1, len(e)):
                        if st[e[y]] == 0:
                            if diff_letters(first[e[x]], first[e[y]]) <= 1:
                                All[e[x]] += 1
                                st[e[y]] += 1
                                removed[e[x]].append(e[y])
                                r.append(e[y])
        return r

    trash = mapf()
    ss = list(set(trash))
    print("\nnumber of kmers to be deleted: {0:d}".format(len(ss)))

    var = pd.DataFrame.from_dict(All, orient='index', columns=['variance'])
    final = kmers1.drop(columns=['tmpid'])
    final['var'] = var['variance']
    final = final.sort_values(by=['var'], ascending=False)
    final = final.drop(trash)
    print("Third step Ends\n")
    final.to_csv(name, index=False)
    print("\nthe final number of kmers: {0:d}".format(len(final)))
