
import tqdm
import networkx as nx
from bs4 import BeautifulSoup
from matplotlib import pyplot as plt
import os
from files_handler import *
from prerequest import *
from distances_and_comparing_functions import *
from matplotlib import pyplot as plt
from spike_analyzing_handler import *
from threading import Thread

def create_amino_acid_netwrok():
    os.chdir(str(directory + "\\") + 'Amino_acid')
    AA_G = nx.Graph(name='AA Network')
    for key in amino_acids_kmers_nets_G.keys():
        AA_G.clear()
        for a_node in amino_acids_kmers_nets_G.get(key):
            for b_node in amino_acids_kmers_nets_G.get(key):
                if calculateAminoAcidsSimilarty(a_node.sub_seq, b_node.sub_seq):
                    AA_G.add_edge(a_node, b_node, weight=calculateAminoAcidsSimilarty(a_node.sub_seq, b_node.sub_seq))

        draw_graph3(AA_G, output_filename=str(key) + '_AA_G_graph_output',
                    show_buttons=True,
                    notebook=False)
        show_Nuclotides_netwrok_stats(AA_G, str(key) + '_NC_G_graph_output')
    plot_spikes(str(directory))
    os.chdir('..')



def create_nuclotides_netwrok():
    print('start the nucoltides node positions')
    nucoltides_kmers_nets_G = {}
    print(len(nucoltides_kmers_nets))
    for key in nucoltides_kmers_nets:
        threads = []
        nucoltides_kmers_nets_G.clear()
        nucoltides_kmers_nets_G = nx.Graph(name='nuclotides kmers_of_position ' + str(key))
        nodes_list = nucoltides_kmers_nets.get(key)
        for node1 in nodes_list:
            thread = Thread(target=compare_nucoltides, args=(node1, nodes_list, nucoltides_kmers_nets_G, key))
            threads.append(thread)
        print("starting with the start")
        for x in tqdm(threads):
            x.start()
        print("starting with the join")
        for x in tqdm(threads):
            x.join()
        print("printing network")
        os.chdir(pajek)
        nx.write_pajek(nucoltides_kmers_nets_G, str(key) + '_NC_G_graph_output1.net')
        os.chdir('..')
        make_netwrok_labels(nucoltides_kmers_nets_G, key)
        os.chdir(nucleotides)
        draw_kmers_graph3(nucoltides_kmers_nets_G,
                          output_filename=str(key) + '_NC_G_graph_output',
                          show_buttons=True,
                          notebook=False)

        show_Nuclotides_netwrok_stats(nucoltides_kmers_nets_G, str(key) + '_NC_G_graph_output')
        os.chdir('..')

def make_netwrok_labels(networkx_graph, key):
    G_labels = nx.Graph(name='AA Netwrok')
    for node, node_attrs in tqdm(networkx_graph.nodes(data=True)):
        G_labels.add_node(node.sub_seq, title=str(str(node.sub_seq) + str(node.subject)),
                          label=str(str(node.sub_seq) + str(node.subject)))
    #         print(node,node_attrs)
    # for each edge and its attributes in the networkx graph
    for source, target, edge_attrs in tqdm(networkx_graph.edges(data=True)):
        # if value/width not specified directly, and weight is specified, set 'value' to 'weight'
        if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
            # place at key 'value' the weight of the edge
            edge_attrs['value'] = edge_attrs['weight']
            G_labels.add_edge(source.sub_seq, target.sub_seq, weight=int(edge_attrs['weight']),
                              title=str(edge_attrs['weight']))
    os.chdir(Amino_acid)
    nx.write_pajek(G_labels, str(key) + '_labels_output1.net')
    os.chdir('..')


def show_Nuclotides_netwrok_stats(networkx_graph, network_name):
    max = 0
    the_key = None
    result = ""
    sequences = ""
    for key in aa_count:
        result += "\n"
        result += (' ' + key + '  ')
        lst = aa_count.get(key)
        result += (str(lst[0]) + '  ' + ' '.join(str(x.sub_seq) for x in lst[1:]))
        if lst[0] > max:
            max = lst[0]
            sequences = aa_count.get(key)
            the_key = key
    result += "\n"
    result += (str(the_key) + '   The max encountered sub sequence is ' + str(max) + '  ' + ' '.join(
        str(x.sub_seq) + ' ' + str(x.position) for x in sequences[1:]))
    most_connected_node_counter = [0, None]
    for node, node_attrs in tqdm(networkx_graph.nodes(data=True)):
        if networkx_graph.degree(node) > most_connected_node_counter[0]:
            most_connected_node_counter[1] = node
            most_connected_node_counter[0] = networkx_graph.degree(node)
    result += "\n"
    result += str(most_connected_node_counter[1]) + ' with this count of connectivity' + str(
        most_connected_node_counter[0])
    with open(str(network_name) + "network_stats.txt", "w") as file:
        file.write(str(result))
    return result

def draw_graph3(networkx_graph, notebook=True, output_filename='empgraph.html', show_buttons=False,
                only_physics_buttons=False):
    from pyvis import network as net
    actual_total_appearances = 0
    real_nodes_count = 0
    real_edges_count = 0
    nuclotide_count_dict = {}
    # make a pyvis network
    pyvis_graph = net.Network(height='500px', width='500px', notebook=notebook)
    pyvis_graph.width = '1000px'
    pyvis_graph = net.Network(height='900px', width='900px', bgcolor='#222222', font_color='white')
    switcher = {
        '1': 'red',
        '2': 'green',
        '3': 'red',
        '4': 'magenta',
        '5': 'cyan',
        '6': 'blue',
        '7': 'black',
        '8': '#eeefff'
    }
    # for each node and its attributes in the networkx graph
    for node, node_attrs in tqdm(networkx_graph.nodes(data=True)):
        color1 = switcher.get(str(node.subject))
        x = nuclotide_count_dict.get(node.sub_seq)
        if x:
            title = "The node Sub sequence is " + str(node.sub_seq)
            if x[0] >= 2:
                x[0] += 1
                nuclotide_count_dict[node.sub_seq] = x
                x.append(node)

                for node1 in x[1:]:
                    title += str(x[0]) + " the subject id is : " + str(
                        node1.subject) + " The sequence id is :" + str(node1.seq_id) + "\n"
                if node.sub_seq not in pyvis_graph.nodes:
                    pyvis_graph.add_node(node.sub_seq, title=title,
                                         label=str(str(node.sub_seq)), color=color1, )
                    real_nodes_count += 1
            else:
                x[0] += 1
                nuclotide_count_dict[node.sub_seq] = x
                x.append(node)
                if x[0] >= 2:
                    for node1 in x[1:]:
                        title += str(x[0]) + " the subject id is : " + str(
                            node1.subject) + " The sequence id is :" + str(node1.seq_id) + "\n"
                    if node.sub_seq not in pyvis_graph.nodes:
                        pyvis_graph.add_node(node.sub_seq, title=title,
                                             label=str(str(node.sub_seq)), color=color1, )
                        real_nodes_count += 1
                        actual_total_appearances += x[0]
        else:
            x = [0]
            x[0] = 1
            x.append(node)
            nuclotide_count_dict[node.sub_seq] = x
    #         print(node,node_attrs)

    # for each edge and its attributes in the networkx graph
    for source, target, edge_attrs in tqdm(networkx_graph.edges(data=True)):
        # if value/width not specified directly, and weight is specified, set 'value' to 'weight'
        if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
            # place at key 'value' the weight of the edge
            edge_attrs['value'] = edge_attrs['weight']
        # add the edge
        x = nuclotide_count_dict.get(source.sub_seq)
        y = nuclotide_count_dict.get(target.sub_seq)
        if x:
            if x[0] >= 2:
                if y:
                    if y[0] >= 2:
                        pyvis_graph.add_edge(source.sub_seq, target.sub_seq, title=str(edge_attrs['weight']),
                                             value=int(edge_attrs['weight']), weight=int(edge_attrs['weight'])
                                             , color='yellow')
                        real_edges_count += 1

    # turn buttons on
    if show_buttons:
        if only_physics_buttons:
            pyvis_graph.show_buttons(filter_=["physics"])
        else:
            pyvis_graph.show_buttons()

    # return and also save
    pyvis_graph.save_graph(output_filename + ".html")
    soup = BeautifulSoup(open(output_filename + ".html"), 'html.parser')
    extra_html = '''
    <div class="record">
        <div class="header">
        <h2>
        '''
    extra_html += "real nodes count =" + str(real_nodes_count) + "\n"
    extra_html += "actual total appearances =" + str(actual_total_appearances) + "\n"
    extra_html += "real edges count =" + str(real_edges_count) + "\n"
    extra_html += "</h2> <div class='title'>"
    extra_html += nx.info(networkx_graph)
    #######################################
    extra_html += "Network density:" + str(nx.density(networkx_graph)) + "\n"
    extra_html += ""
    # plt.plot(list(nuclotide_count_dict.keys()), [nuclotide_count_dict.get(x)[0] for x in nuclotide_count_dict.keys()])
    # plt.savefig("tempfig.png")
    # plt.show()

    extra_html += '''
            </div>
        </div>
    </div>'''

    div = soup.select("#config")
    soup.find_all('div', {"id": "mynetwork"})[-1].insert_after(BeautifulSoup(extra_html, 'html.parser'))
    with open(output_filename + ".html", "w") as file:
        file.write(str(soup))


def draw_kmers_graph3(networkx_graph, notebook=True, output_filename='empgraph', show_buttons=False,
                      only_physics_buttons=False):
    # import
    from pyvis import network as net
    real_nodes_count = 0
    real_edges_count = 0
    # make a pyvis network
    pyvis_graph = net.Network(height='900px', width='900px', bgcolor='#222222', font_color='white')
    switcher = {
        '1': 'red',
        '2': 'green',
        '3': 'red',
        '4': 'magenta',
        '5': 'cyan',
        '6': 'blue',
        '7': 'black',
        '8': '#eeefff'
    }
    nuclotide_count_dict = {}
    # for each node and its attributes in the networkx graph
    nodes_count = len(networkx_graph.nodes(data=True))
    for node, node_attrs in tqdm(networkx_graph.nodes(data=True)):
        color1 = switcher.get(str(node.subject))
        x = nuclotide_count_dict.get(node.sub_seq)
        if x:

            title = "The node Sub sequence is " + str(node.sub_seq)
            if x[0] >= 2:
                x[0] += 1
                nuclotide_count_dict[node.sub_seq] = x
                x.append(node)

                for node1 in x[1:]:
                    title += " the subject id is : " + str(
                        node1.subject) + " The sequence id is :" + str(node1.seq_id) + "\n"
                pyvis_graph.add_node(node.sub_seq, title=title, label=str(str(node.sub_seq)), color=color1, )
            else:
                x[0] += 1
                nuclotide_count_dict[node.sub_seq] = x
                x.append(node)
                if x[0] >= 2:
                    for node1 in x[1:]:
                        title += " the subject id is : " + str(
                            node1.subject) + " The sequence id is :" + str(node1.seq_id) + "\n"
                    if node.sub_seq not in pyvis_graph.nodes:
                        pyvis_graph.add_node(node.sub_seq, title=title,
                                             label=str(str(node.sub_seq)), color=color1, )
                        real_nodes_count += 1

        else:
            x = [0]
            x[0] = 1
            x.append(node)
            nuclotide_count_dict[node.sub_seq] = x

    #         print(node,node_attrs)

    # for each edge and its attributes in the networkx graph
    for source, target, edge_attrs in tqdm(networkx_graph.edges(data=True)):

        # if value/width not specified directly, and weight is specified, set 'value' to 'weight'
        if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
            # place at key 'value' the weight of the edge
            edge_attrs['value'] = edge_attrs['weight']
        # add the edge
        title = ""
        x = nuclotide_count_dict.get(source.sub_seq)
        y = nuclotide_count_dict.get(target.sub_seq)
        if x:
            if x[0] >= 2:
                if y:
                    if y[0] >= 2:
                        pyvis_graph.add_edge(source.sub_seq, target.sub_seq, title=str(3 - int(edge_attrs['weight'])),
                                             value=3 - int(edge_attrs['weight']), weight=3 - int(edge_attrs['weight']))
                        real_edges_count += 1

        # turn buttons on
    if show_buttons:
        if only_physics_buttons:
            pyvis_graph.show_buttons(filter_=["physics"])
        else:
            pyvis_graph.show_buttons()
    # nx.number_weakly_connected_components(G)
    if nx.number_connected_components(networkx_graph):
        nx.number_connected_components(networkx_graph)

    pyvis_graph.save_graph(output_filename + ".html")
    soup = BeautifulSoup(open(output_filename + ".html"), 'html.parser')
    extra_html = '''
       <div class="record">
           <div class="header">
           <h2>
           '''
    extra_html += "real nodes count =" + str(real_nodes_count) + "\n"
    extra_html += "real edges count =" + str(real_edges_count) + "\n"
    extra_html += "</h2> <div class='title'>"
    extra_html += nx.info(networkx_graph)
    #######################################
    extra_html += "Network density:" + str(nx.density(networkx_graph)) + "\n"
    extra_html += "" + plt.plot(nuclotide_count_dict.keys(), nuclotide_count_dict.values)

    extra_html += '''
               </div>
           </div>
       </div>'''
    div = soup.select("#config")
    soup.find_all('div', {"id": "mynetwork"})[-1].insert_after(BeautifulSoup(extra_html, 'html.parser'))
    with open(output_filename + ".html", "w") as file:
        file.write(str(soup))
