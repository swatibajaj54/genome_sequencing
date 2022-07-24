import itertools
import math


def string_composition(text, k):
    """returns string of k length from the input text in any sequence
    (eg: text- atgg, k =2, output = at, tg, gg"""
    i = 0
    string_comp = []
    while i <= len(text) - k:
        string = text[i: i + k]
        string_comp.append(string)
        i = i + 1
    return string_comp


def string_spelled_by_genome_path(genome_path):
    """genome path is string composition in sequence.
    The function returns string_text from the genome path"""
    k = len(genome_path[0])
    reconstruction = [genome_path[0][0: k - 1]]
    for item in genome_path:
        reconstruction.append(item[k - 1])
    string_reconstruction = "".join(reconstruction)
    return string_reconstruction


def overlap_graph_problem(patterns):
    """Input: A collection Patterns of k-mers.
    Output: The overlap graph Overlap(Patterns), in the form of an adjacency list. (in any order)"""
    k = len(patterns[0])
    prefix_dict = {}
    adjacency_list = {}
    for item in patterns:
        key = item[0: k - 1]
        val = prefix_dict.get(key)
        val = val + " " if val else ""
        val = val + item
        prefix_dict[key] = val
    for key in patterns:
        suffix = key[1:k]
        if prefix_dict.get(suffix):
            adjacency_list[key] = prefix_dict.get(suffix)
    return adjacency_list


def de_bruijn_graph(patterns):
    """returns DeBruijnk(Text), in the form of an adjacency list."""
    de_bruijn_dict = {}
    k = len(patterns[0])
    for item in patterns:
        key = item[0: k - 1]
        value = item[1:k]
        val = de_bruijn_dict.get(key)
        if val:
            val += " "
        else:
            val = ""
        val = val + value
        de_bruijn_dict[key] = val
    for key in de_bruijn_dict:
        val = de_bruijn_dict[key].split()
        de_bruijn_dict[key] = val
    return de_bruijn_dict


def de_bruijn_graph_from_text(text, k):
    """Input: An integer k and a string Text.
    Output: DeBruijnk(Text), in the form of an adjacency list."""
    patterns = string_composition(text, k)
    de_bruijn_dict = de_bruijn_graph(patterns)
    return de_bruijn_dict


def eulerian_cycle_3(adjacency_dict: dict):
    """returns eulerian cycle (close circuit) (each edge visit only once) from adjacency_list"""
    edges_list = get_edges_list_from(adjacency_dict)
    stack = []
    res = []
    keys = [*adjacency_dict]

    stack.append(keys[0])
    while len(stack) != 0:
        top_elem = stack[-1]
        vertices = adjacency_dict[top_elem]

        # process each vertex
        for vertex in vertices:
            edge = f"{top_elem}-{vertex}"
            if edges_list[edge] == "not_visited":
                stack.append(vertex)
                edges_list[edge] = "visited"
                break
        else:
            processed_node = stack.pop()
            res.append(processed_node)

    res.reverse()
    return res


def get_edges_list_from(adjacency_dict: dict):
    """returns dictionary of edges e.g: path 0:1 (0-1 = not_visited)"""
    res: dict = {}
    for key in adjacency_dict:
        arr = adjacency_dict[key]
        for vertex in arr:
            res[f"{key}-{vertex}"] = "not_visited"
    return res


def incoming_vertices(adjacency_dict):
    """returns a dictionary which is reverse of adjacency_dict with respect to its key value pair"""
    incoming_vertices_dict = {}
    for key in adjacency_dict:
        for val in adjacency_dict[key]:
            incoming_vertices_dict[val] = incoming_vertices_dict.get(val, []) + [key]
    return incoming_vertices_dict


def eulerian_path(adjacency_dict: dict):
    """returns eulerian path (each edge visit only once) from adjacency_list"""
    stack = []
    res = []
    in_vertices = incoming_vertices(adjacency_dict)
    odd_vertex = ''
    for key in adjacency_dict:
        if key not in in_vertices:
            odd_vertex = key
        else:
            outgo_vertices = adjacency_dict[key]
            for vertex in in_vertices:
                if vertex == key:
                    incom_vertices = in_vertices[vertex]
                    if len(outgo_vertices) != len(incom_vertices):
                        odd_vertex = key
    stack.append(odd_vertex)
    new_edge = ''
    for key in in_vertices:
        if key not in adjacency_dict:
            new_edge = key
    adjacency_dict[new_edge] = [odd_vertex]
    # print(adjacency_dict)
    edges_list = get_edges_list_from(adjacency_dict)
    edges_list[f"{new_edge}-{odd_vertex}"] = "visited"
    # print("Edges list: ", edges_list)
    while len(stack) != 0:
        top_elem = stack[-1]
        vertices = adjacency_dict[top_elem]

        # process each vertex
        for vertex in vertices:
            edge = f"{top_elem}-{vertex}"
            if edges_list[edge] == "not_visited":
                stack.append(vertex)
                edges_list[edge] = "visited"
                break
        else:
            processed_node = stack.pop()
            # print(processed_node)
            res.append(processed_node)

    res.reverse()
    return res


def string_reconstruction(patterns):
    """returns text from patterns(a list) which are puzzled """
    adjacency_dict = de_bruijn_graph(patterns)
    # print('adjacency_list:', adjacency_dict)
    for key in adjacency_dict:
        val = adjacency_dict[key].split()
        adjacency_dict[key] = val
    path = eulerian_path(adjacency_dict)
    # print('path:', path)
    text = string_spelled_by_genome_path(path)
    # print('text:', text)
    return text


def k_universal_circular_string(k):
    """returns universal binary string of length k"""
    binary_kmers = set()
    for seq in itertools.product("01", repeat=k):
        binary_kmers.add("".join(seq))
    # print(binary_kmers)
    adjacency_dict = de_bruijn_graph(list(binary_kmers))
    # print('dict:', adjacency_dict)
    circuit = eulerian_cycle_3(adjacency_dict)
    # print(circuit)
    universal_string = string_spelled_by_genome_path(circuit)
    # print(universal_string)
    valid_length = 2 ** k
    circular_universal_string = universal_string[:valid_length]
    return circular_universal_string


def read_pair_composition(text, k, d):
    """returns read pairs of length k each from text with a gap length d """
    read_pairs = []
    i = 0
    while i <= len(text) - (2 * k + d):
        first_kmer = text[i:i + k]
        second_kmer = text[i + k + d:i + k + d + k]
        read_pair = (first_kmer, "|", second_kmer)
        read_pairs.append(read_pair)
        i = i + 1
    read_pairs.sort()
    return read_pairs


def string_reconstruction_from_paired_de_bruijn_graph_path(text, d, k):
    """contruct a genome with paired de bruijn graph path"""
    prefix_pattern = []
    suffix_pattern = []
    for item in text:
        prefix_pattern.append(item[0])
        suffix_pattern.append(item[2])
    prefix_string = string_spelled_by_genome_path(prefix_pattern)
    # print(prefix_string)
    suffix_string = string_spelled_by_genome_path(suffix_pattern)
    # print(suffix_string)
    i = k + d
    item = suffix_string[i - k - d]
    item2 = prefix_string[i]
    # print(item2)
    # print(item)
    while i <= len(prefix_string) - k - d:
        if prefix_string[i] != suffix_string[i - k - d]:
            return "there is no string spelled by the gapped patterns"
        i = i + 1
    return prefix_string + suffix_string[len(suffix_string) - k - d:]


def sring_reconstruction_from_paired_reads(patterns, d):
    """returns genome string from read pairs."""
    de_bruijn_dict = {}
    k = len(patterns[0][0])
    for item in patterns:
        key = item[0][: k - 1], "|", item[2][: k - 1]
        value = [(item[0][1:k], "|", item[2][1:k])]
        val = de_bruijn_dict.get(key)
        if val:
            val = value.append(value)
        else:
            val = value
        de_bruijn_dict[key] = val
    # print(de_bruijn_dict)
    path = eulerian_path(de_bruijn_dict)
    # print(path)
    string = string_reconstruction_from_paired_de_bruijn_graph_path(path, d, k)
    return string


def max_non_branching_path(adjacency_dict):
    paths = []
    in_vertices = incoming_vertices(adjacency_dict)
    for key in adjacency_dict:
        if key not in in_vertices or len(adjacency_dict[key]) > 1 or len(in_vertices[key]) > 1:
            if len(adjacency_dict[key]) > 0:
                for out_edge in adjacency_dict[key]:
                    non_branching_path = []
                    non_branching_path.extend([key, out_edge])
                    if out_edge in adjacency_dict:
                        if len(adjacency_dict[out_edge]) == 1 and out_edge in in_vertices:
                            if len(in_vertices[out_edge]) == 1:
                                non_branching_path.extend(adjacency_dict[
                                                              out_edge])  # extend will not work with dna string find alternative to convert list
                    paths.append(non_branching_path)

    for key in adjacency_dict:
        cycle = []
        used_key = ''
        if key in in_vertices and adjacency_dict[key] == in_vertices[key]:
            if key == used_key:
                cycle.extend(adjacency_dict[key])
                cycle.extend([key])
                cycle.extend(adjacency_dict[key])
                used_key = key
                paths.append(cycle)
                # if key == adjacency_dict[used_key]:
                # continue
    print(paths)
