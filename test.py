## Initialize PyIR and set example file for processing
import gzip
import json
#
# from crowelab_pyir import PyIR
#
# FILE = './covid_vaccine_newallSelect Start And End Position11225510000cdr3_fasta'+'.fasta'
#
# pyirfile = PyIR(query=FILE,)
# result = pyirfile.run()
#
# # Prints the output file
# print(result)
#
# # with gzip.open(str('covid_vaccine_newallSelect Start And End Position11225510000cdr3_fasta')+".json.gz", "r") as f:
# #     data = f.read()
# #     j = json.loads(data.decode('utf-8'))
# #     print(type(j))
## Initialize PyIR and set example file for processing
import pandas as pd
from crowelab_pyir import PyIR
from Bio import SeqIO
FILE = './covid_vaccine_new3Select Start And End Position11215610000cdr3_fasta.fasta'

pyirfiltered = PyIR(query=FILE, args=['--outfmt', 'dict'])
result = pyirfiltered.run()
df = pd.DataFrame.from_dict(result, orient = 'index')

#Prints size of Python returned dictionary
print(df)
with open("api_results.csv", 'w', newline='') as new_file:
    df.to_csv(path_or_buf=new_file)