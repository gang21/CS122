import csv
from typing import List, Dict, Iterable
import random
from collections import defaultdict, Counter
from copy import deepcopy
import random
from venv import create

def read_spectrum(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()

    reads = []
    for read in lines:
        if read[0] != ">":
            read = read.strip()
            reads.append(read)
    return reads

def breakup(reads):
    split_dict = {}
    for read in reads:
        first = read[:-1]
        second = read[1:]
        split_dict[first] = second
    return split_dict

def find_start(split_dict):
    potentials = []
    for key,value in split_dict.items():
        if key not in split_dict.values() and value in split_dict.keys():
            potentials.append(key)
    print(len(potentials))
    return potentials

def recreate_path(reads, splits, start):
    start_kmer = start + splits[start][-1]
    val = splits[start]
    path = [start_kmer]
    index = reads.index(start_kmer)
    indexes = [index]

    while val in splits.keys():
        kmer = splits[val]
        full_kmer = val + kmer[-1]
        path.append(full_kmer)
        index = reads.index(full_kmer)
        indexes.append(index)
        # print(kmer, "  -  ", index)
        val = kmer
    positions = []
    for position in indexes:
        final = ">read_" + str(position)
        positions.append(final)
    return positions

def create_csv(values):
      with open("predictions.csv", 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        for value in values:
            csv_writer.writerow([value])

def genome_path(path):
    """Forms the genome path formed by a collection of patterns."""
    final_string = path[0]
    for i in range(1,len(path)):
        read = path[i]
        final_string = final_string + read[-1]
    return final_string

def de_bruijn_kmers(k_mers):
    """Forms the de Bruijn graph of a collection of k-mers."""
    final_dict = {}
    for i in range(len(k_mers)):
        key = k_mers[i][:-1]
        vals = []
        for j in range(len(k_mers)):
            if k_mers[j][:-1] == key:
                vals.append(k_mers[j][1:])
        final_dict[key] = vals
    return final_dict

def find_cycle(graph, start):
    cycle = []
    u = graph[start].pop()
    while u != start:
        cycle.append(u)
        u = graph[u].pop()
    cycle.append(u)

    # Clean up nodes which have no edges.
    toRemove = [k for k, v in graph.items() if not v]
    for k in toRemove:
        del graph[k]

    return cycle

def find_eulerian_cycle(graph, start=0):
    cycle = [start] + find_cycle(graph, start)
    updated = True
    while updated:
        updated = False
        for i, start in enumerate(cycle):
            # If an edge starting from the node still exists,
            # insert new cycle.
            if start in graph:
                updated = True
                cycle = cycle[:i+1] + find_cycle(graph, start) + cycle[i+1:]
                break

    return cycle

def add_imaginary_edge(graph):
    outgoingEdgeCounts, incomingEdgeCounts = Counter(), Counter()
    for u in graph:
        outgoingEdgeCounts[u] += len(graph[u])
        for v in graph[u]:
            incomingEdgeCounts[v] += 1 

    start = list((incomingEdgeCounts - outgoingEdgeCounts).keys())[0]
    end = list((outgoingEdgeCounts - incomingEdgeCounts).keys())[0]

    # Add imaginary edge.
    graph[start].append(end)
    return start, end

def find_eulerian_path(graph):
    start, end = add_imaginary_edge(graph)
    cycle = find_eulerian_cycle(graph, start=end)[:-1]
    for i in range(len(cycle)):
        if cycle[i] == start and cycle[(i+1) % len(cycle)] == end:
            path = cycle[i+1:] + cycle[:i+1]

    return path

def eulerian_path(g: Dict[int, List[int]]) -> Iterable[int]:
    start_node = None
    end_node = None
    in_degrees = {}
    out_degrees = {}
    final_path = []
    graph = deepcopy(g)
    for key, val in g.items():
        for items in val:
            if key in out_degrees.keys():
                out_degrees[key] = out_degrees[key] + 1
            else: 
                out_degrees[key] = 1
            if items in in_degrees.keys():
                in_degrees[items] = in_degrees[items] + 1 
            else:
                in_degrees[items] = 1
            if items not in graph.keys():
                graph[items] = []
    # find start and end nodes
    for key, val in in_degrees.items():
        if key not in out_degrees.keys():
            end_node = key
            continue
        in_ = in_degrees[key]
        out_ = out_degrees[key]
        
        if in_ == out_ - 1:
            start_node = key
            break
            
    result = find_eulerian_path(graph)   
    return result


def string_reconstruction(patterns, k):
    dB = de_bruijn_kmers(patterns)
    path = eulerian_path(dB)
    text = genome_path(path)
    
    return text

def main():
    reads = read_spectrum("project2a_spectrum.fasta")
    splits = breakup(reads)
    length = len(reads[0])
    reconstructed_genome = string_reconstruction(reads, length)
    print(reconstructed_genome)
    positions = []
    for i in range(len(reconstructed_genome)-length):
        kmer = reconstructed_genome[i:i+length]
        index = reads.index(kmer)
        full_read = ">read_" + str(index)
        positions.append(full_read)
    create_csv(positions)
    
    
if __name__ == "__main__":
    main()
