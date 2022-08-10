import networkx as nx

global kmers_nets
kmers_nets = {}
global nucoltides_kmers_nets
nucoltides_kmers_nets = {}
nucoltides_kmers_nets_edges_count = {}
replaced_nucoltides_kmers_nets = {}
replaced_nucoltides_kmers_nets_edges_count = {}
amino_acid_kmers_nets = {}

G = nx.Graph(name='CDR3 from DB')
nucoltides_kmers_nets_G = {}
amino_acids_kmers_nets_G = {}

NC_G = nx.Graph(name='AA Netwrok')
aa_count = {"TTT": [0], "CTT": [0], "ATT": [0], "GTT": [0],
            "TTC": [0], "CTC": [0], "ATC": [0], "GTC": [0],
            "TTA": [0], "CTA": [0], "ATA": [0], "GTA": [0],
            "TTG": [0], "CTG": [0], "ATG": [0], "GTG": [0],
            "TCT": [0], "CCT": [0], "ACT": [0], "GCT": [0],
            "TCC": [0], "CCC": [0], "ACC": [0], "GCC": [0],
            "TCA": [0], "CCA": [0], "ACA": [0], "GCA": [0],
            "TCG": [0], "CCG": [0], "ACG": [0], "GCG": [0],
            "TAT": [0], "CAT": [0], "AAT": [0], "GAT": [0],
            "TAC": [0], "CAC": [0], "AAC": [0], "GAC": [0],
            "TAA": [0], "CAA": [0], "AAA": [0], "GAA": [0],
            "TAG": [0], "CAG": [0], "AAG": [0], "GAG": [0],
            "TGT": [0], "CGT": [0], "AGT": [0], "GGT": [0],
            "TGC": [0], "CGC": [0], "AGC": [0], "GGC": [0],
            "TGA": [0], "CGA": [0], "AGA": [0], "GGA": [0],
            "TGG": [0], "CGG": [0], "AGG": [0], "GGG": [0],
            "---": [0]
            }
aa_most_frequent = {}
spike_dic = {}
amino_acid_list = ['ATA', 'ATC', 'ATT', 'ATG',
                   'ACA', 'ACC', 'ACG', 'ACT',
                   'AAC', 'AAT', 'AAA', 'AAG',
                   'AGC', 'AGT', 'AGA', 'AGG',
                   'CTA', 'CTC', 'CTG', 'CTT',
                   'CCA', 'CCC', 'CCG', 'CCT',
                   'CAC', 'CAT', 'CAA', 'CAG',
                   'CGA', 'CGC', 'CGG', 'CGT',
                   'GTA', 'GTC', 'GTG', 'GTT',
                   'GCA', 'GCC', 'GCG', 'GCT',
                   'GAC', 'GAT', 'GAA', 'GAG',
                   'GGA', 'GGC', 'GGG', 'GGT',
                   'TCA', 'TCC', 'TCG', 'TCT',
                   'TTC', 'TTT', 'TTA', 'TTG',
                   'TAC', 'TAT', 'TAA', 'TAG',
                   'TGC', 'TGT', 'TGA', 'TGG', ]

proteinNumber = 64
countrow = 0