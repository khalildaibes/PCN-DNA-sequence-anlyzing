import pandas as pd
from tqdm import tqdm

tqdm.pandas()
from collections import Counter


# ID=10
# path=r"/home/saif/metadata/covid/covid{}_seqK.csv".format(ID)
# name="/home/saif/zak/cov1/cov{}_byScore.csv".format(ID)
# col_list = ["k-mer","ai","clone_id"]
# kmerss=pd.read_csv(path,usecols=col_list)

def uniqueK(kmers, name):
    print("Second step (extract the unique k-mers):\n")

    # print("unique starts here")
    kmers = kmers.rename(columns={'Unique-SeqKmer': 'Kmers'}, inplace=False)
    kmers = kmers.rename(columns={'k-mer': 'Kmers'}, inplace=False)

    print("number of k-mers before:" + str(len(kmers)))
    kmers['id'] = kmers.index

    l = kmers.values.tolist()

    d = {}
    for i in tqdm(l):
        d[i[0]] = []
    for j in tqdm(l):
        if (j[3] not in d[j[0]]):
            d[j[0]].append(j[3])

    new_df = pd.DataFrame(kmers['Kmers'])
    l = new_df.values.tolist()
    flat_list = []
    for sublist in l:
        for item in sublist:
            flat_list.append(item)
    l = flat_list

    def count_uniqe(lst):
        new_vals = Counter(l).most_common()
        new_vals = new_vals[::1]  # this sorts the list in scending order
        return new_vals

    new_list = count_uniqe(l)

    df = pd.DataFrame(new_list, columns=['kmer', 'score'])

    print("number of unique k-mers: " + str(len(df)))

    def inx(kmer):
        return d[kmer][0]

    df['id'] = df.kmer.apply(inx)
    df.to_csv(name, index=False)
    print("Second step Ends\n")





