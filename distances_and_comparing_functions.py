
def compare_nucoltides(node1, nodes_list, nucoltides_kmers_nets_G, key):
    for node2 in nodes_list:
        if node1.sub_seq != node2.sub_seq and node1.sub_seq[0:2] in node2.sub_seq or node1.sub_seq[3:5] in \
                node2.sub_seq or node1.sub_seq[6:8] in node2.sub_seq:
            w = calculateHammingDistance(node1.sub_seq, node2.sub_seq)
            if w:
                nucoltides_kmers_nets_G.add_edge(node1, node2, weight=w)
def calculateHammingDistance(protein1, protein2):
    count = 0
    rangeofseq = min(len(protein1), len(protein2))
    for j in range(0, rangeofseq):
        if protein1[j] == protein2[j]:
            pass
        else:
            count = count + 1
    if count > 2:
        return 0

    return count


def calculateAminoAcidsSimilarty(protein1, protein2):
    count = 0
    rangeofseq = min(len(protein1), len(protein2))
    for j in range(0, rangeofseq):
        if protein1[j] == protein2[j]:
            pass
        else:
            count = count + 1
    if count > 1:
        return 0
    return count