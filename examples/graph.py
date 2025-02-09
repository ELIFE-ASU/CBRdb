import os

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

import CBRdb


def get_compound_list(path_c="../data/kegg_data_C.csv", head=None):
    """
    Reads the compound data from a CSV file and returns a list of compound IDs.

    Parameters:
    path_c (str): The file path to the compound data CSV file. Defaults to "../data/kegg_data_C.csv".
    head (int, optional): The number of rows to read from the CSV file. If None, all rows are read. Defaults to None.

    Returns:
    list: A list of compound IDs from the CSV file.
    """
    data_c = pd.read_csv(os.path.abspath(path_c))
    if head is not None:
        data_c = data_c.head(head)
    return data_c["compound_id"].tolist()


def generate_nodes_from_compounds(compound_list, graph=None):
    """
    Generates nodes in the graph based on the given list of compounds.

    Parameters:
    compound_list (list): A list of compound IDs.
    graph (nx.Graph, optional): A NetworkX graph where nodes will be added. If None, a new graph is created. Defaults to None.

    Returns:
    nx.Graph: The updated graph with nodes added based on the compounds.
    """
    if graph is None:
        graph = nx.Graph()
    for compound_id in compound_list:
        graph.add_node(compound_id)
    return graph


def get_reaction_list(path_r="../data/kegg_data_R.csv", head=None):
    """
    Reads the reaction data from a CSV file and returns a list of reactions.

    Parameters:
    path_r (str): The file path to the reaction data CSV file. Defaults to "../data/kegg_data_R.csv".
    head (int, optional): The number of rows to read from the CSV file. If None, all rows are read. Defaults to None.

    Returns:
    list: A list of reaction strings from the CSV file.
    """
    data_r = pd.read_csv(os.path.abspath(path_r))
    if head is not None:
        data_r = data_r.head(head)
    return data_r["reaction"].tolist()


def generate_edges_from_reactions(reaction_list, graph, f_print=True):
    """
    Generates edges in the graph based on the given list of reactions.

    Parameters:
    reaction_list (list): A list of reaction strings.
    graph (nx.Graph): A NetworkX graph where edges will be added.
    f_print (bool): A flag to enable or disable printing of debug information.

    Returns:
    nx.Graph: The updated graph with edges added based on the reactions.
    """
    for reaction in reaction_list:
        # Convert the reaction string to reactants and products
        reactants, products = CBRdb.eq_to_dict(reaction)
        if f_print:
            print(f"Reaction: {reaction}", flush=True)
            print(f"Reactants: {reactants}", flush=True)
            print(f"Products: {products}", flush=True)

        # Iterate over each reactant and product to add edges
        for reactant in reactants:
            for product in products:
                # Check if both reactant and product are in the graph
                if reactant in graph.nodes and product in graph.nodes:
                    graph.add_edge(reactant, product)
                    if f_print:
                        print(f"Added edge between {reactant} and {product}", flush=True)
                else:
                    if f_print:
                        print(f"Node {reactant} or {product} not in graph", flush=True)
                    pass
    return graph


def draw_graph(graph):
    # pos = nx.spring_layout(graph)
    # nx.draw(graph, pos, with_labels=True, node_size=700)
    nx.draw(graph, node_size=10)
    plt.show()


def count_disjointed_graphs(graph):
    """
    Counts the number of disjointed (connected) components in the given graph.

    Parameters:
    graph (nx.Graph): A NetworkX graph.

    Returns:
    int: The number of disjointed components in the graph.
    """
    return nx.number_connected_components(graph)


def count_nodes_and_edges(graph):
    """
    Counts the number of nodes and edges in the given graph.

    Parameters:
    graph (nx.Graph): A NetworkX graph.

    Returns:
    tuple: A tuple containing the number of nodes and edges in the graph.
    """
    num_nodes = graph.number_of_nodes()
    num_edges = graph.number_of_edges()
    return num_nodes, num_edges


def get_graph_info(graph):
    """
    Retrieves and prints information about the given graph, including the number of nodes, edges, and disjointed components.

    Parameters:
    graph (nx.Graph): A NetworkX graph.

    Returns:
    tuple: A tuple containing the number of nodes, edges, and disjointed components in the graph.
    """
    num_d = count_disjointed_graphs(graph)
    num_n, num_e = count_nodes_and_edges(graph)
    print(f"N nodes: {num_n}, N edges: {num_e}, N disjointed: {num_d}", flush=True)
    return num_n, num_e, num_d


def prune_nodes_with_no_edges(graph):
    """
    Prunes off nodes that have no edges in the given graph.

    Parameters:
    graph (nx.Graph): A NetworkX graph.

    Returns:
    nx.Graph: The pruned graph with nodes that have no edges removed.
    """
    nodes_to_remove = [node for node in graph.nodes if graph.degree(node) == 0]
    graph.remove_nodes_from(nodes_to_remove)
    return graph


if __name__ == "__main__":
    print("Program started", flush=True)
    compound_list = get_compound_list(head=10)
    graph = generate_nodes_from_compounds(compound_list)

    reaction_list = get_reaction_list(head=10)
    graph = generate_edges_from_reactions(reaction_list, graph)

    draw_graph(graph)
    pruned_graph = prune_nodes_with_no_edges(graph)
    draw_graph(pruned_graph)

    print("Program finished", flush=True)
