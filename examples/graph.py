import networkx as nx
import matplotlib.pyplot as plt

import networkx as nx

def get_compund_list(data):
    return data["compound_id"].tolist()

def generate_nodes_from_compounds(compound_list, graph=None):
    if graph is None:
        graph = nx.Graph()
    for compound_id in compound_list:
        graph.add_node(compound_id)
    return graph




# Example usage
compound_list = ['C00001', 'C00002', 'C00003']
graph = generate_nodes_from_compounds(compound_list)
print(graph.nodes)  # Output: ['C00001', 'C00002', 'C00003']


def draw_graph(graph):
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos, with_labels=True, node_size=700)
    plt.show()

if __name__ == "__main__":
    data_c = pd.read_csv("../data/kegg_data_C.csv.zip")

    G = nx.Graph()
    G.add_node(1)
    G.add_node(2)
    G.add_node(3)
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    G.add_edge(1, 3)
    draw_graph(G)
    print("Program finished", flush=True)