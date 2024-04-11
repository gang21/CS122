from collections import defaultdict, Counter
import math
import csv

def create_csv(values):
    with open("predictions.csv", 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        for value in values:
            csv_writer.writerow([value])

def read_reads(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()

    reads = []
    for read in lines:
        if read[0] != ">":
            read = read.strip()
            reads.append(read)
    return reads

def break_into_kmers(reads, k):
    kmers = []

    for read in reads:
        for i in range(0, len(read) - k + 1):
            kmers.append(read[i:(i + k)])

    return kmers

def get_kmer_frequencies(kmers):
    kmer_to_frequency = {}

    for kmer in kmers:
        if kmer in kmer_to_frequency:
            kmer_to_frequency[kmer] += 1
        else:
            kmer_to_frequency[kmer] = 1

    return kmer_to_frequency

def form_de_bruijn_graph(kmer_to_frequency):
    adjacency_list = {}
    for kmer in kmer_to_frequency:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if prefix in adjacency_list:
            adjacency_list[prefix].append(suffix)
        else:
            adjacency_list[prefix] = [suffix]

    return adjacency_list

def get_node_degrees(adjacency_list):
    node_degrees = {}
    for outgoing_node, incoming_nodes in adjacency_list.items():
        if outgoing_node in node_degrees:
            node_degrees[outgoing_node][1] += len(incoming_nodes)
        else:
            node_degrees[outgoing_node] = [0, len(incoming_nodes)]

        for incoming_node in incoming_nodes:
            if incoming_node in node_degrees:
                node_degrees[incoming_node][0] += 1
            else:
                node_degrees[incoming_node] = [1, 0]

    return node_degrees

def find_paths(adjacency_list):
    node_degrees = get_node_degrees(adjacency_list)
    paths = []

    for v in list(node_degrees.keys()):
        if node_degrees[v] != [1,1]:
            if node_degrees[v][1] > 0:
                for w in adjacency_list[v]:
                    non_branching_path = [v,w]
                    while node_degrees[w] == [1,1]:
                        u = adjacency_list[w][0]
                        non_branching_path.append(u)
                        w = u
                    paths.append(non_branching_path)
    for path in paths:
        for node in path:
            if node in adjacency_list:
                del adjacency_list[node]

    for node in list(adjacency_list.keys()):
        if node_degrees[node] != [1,1]:
            del adjacency_list[node]

    while adjacency_list:
        start_node = list(adjacency_list.keys())[0]
        
        curr_node = start_node
        next_node = adjacency_list[start_node][0]

        cycle = [start_node]

        first_time_visiting_start_node = True
        while curr_node != start_node or first_time_visiting_start_node:
            first_time_visiting_start_node = False

            del adjacency_list[curr_node]
            cycle.append(next_node)

            curr_node = next_node

            if next_node not in adjacency_list:
                continue
            next_node = adjacency_list[next_node][0]
        
        paths.append(cycle)

    return paths

def form_contigs(paths):
    contigs = []

    for path in paths:
        contig = path[0]
        for i in range(1, len(path)):
            contig += path[i][-1]
        contigs.append(contig)

    return contigs

def smith_waterman(v, w):
    n = len(v)
    m = len(w)
    s = [[0 for i in range(m+1)] for j in range(n+1)]
    Backtrack = [[0 for i in range(m+1)] for j in range(n+1)]
    for i in range(n+1):
        s[i][0] = -1 * i
    for j in range(m+1):
        s[0][j] = -1 * j
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = 0
            if v[i-1] == w[j-1]:
                match = 3
            val = max(s[i-1][j]-1, s[i][j-1]-1, s[i-1][j-1]+match)
            s[i][j] = val
    return s[i][j]

def create_hash_table(reference_genome, kmer_length):
    n = len(reference_genome) - kmer_length + 1
    hash_table = defaultdict(list)
    for i in range(n):
        kmer = reference_genome[i:i+kmer_length]
        hash_table[kmer].append(i)
    return hash_table

if __name__ == "__main__":
    input_reads = read_reads("project2b_reads.fasta")
    kmers = break_into_kmers(input_reads, 15)
    kmer_to_frequency = get_kmer_frequencies(kmers)
    kmer_to_frequency = {kmer: freq for kmer, freq in kmer_to_frequency.items() if freq > 3}
    adjacency_list = form_de_bruijn_graph(kmer_to_frequency)
    paths = find_paths(adjacency_list)
    contigs = form_contigs(paths)

    longest = ''.join(contigs)
    # longest = sorted(contigs, key=len, reverse=True)
    # longest = longest[0]
    # # longest = longests[0]
    print(len(longest))
    maxes_per_read = {}
    hash_table = create_hash_table(longest, 16)
    for i in range(len(input_reads)):
        read = input_reads[i]
        len_read = len(read)
        potential_kmers = []
        val1 = read[0:16]
        val2 = read[16:32]
        val3 = read[32:48]
        if val1 in hash_table.keys():
            potential_kmers = potential_kmers + hash_table[val1]
        if val2 in hash_table.keys():
            potential_kmers = potential_kmers + [x - 16 for x in hash_table[val2]]
        if val3 in hash_table.keys():
            potential_kmers = potential_kmers + [x - 32 for x in hash_table[val3]]
        maxes = []
        for position in potential_kmers:
            genome = longest[position:position+len_read]
            sw = smith_waterman(genome, read)
            maxes.append({sw:position})
        if len(maxes) != 0:
            try:
                maximum = max(maxes)
            except:
                maximum = maxes[0]
            # print("Genome: ", i, "  |  Read: ", math.floor(j/2), " -- ", list(maximum.values())[0])
            key = read
            max_val = list(maximum.keys())[0]
            max_pos = list(maximum.values())[0]
            if key in maxes_per_read.keys():
                read_dict = maxes_per_read[key]
                get_max = read_dict['max']
                if max_val > get_max:
                    value = {"max": max_val, "position": max_pos}
                    maxes_per_read[key] = value
            else:
                # print("we r in the else")
                value = {"max": max_val, "position": max_pos}
                maxes_per_read[key] = value

    # sort dictionary and get reads
    sorted_maxes_per_read = dict(sorted(maxes_per_read.items(), key=lambda item: item[1]['position']))
    # print(sorted_maxes_per_read)

    final_list = []
    for key in sorted_maxes_per_read.keys():
        if sorted_maxes_per_read[key]['max'] > 100:
            entry = ">read_" + str(input_reads.index(key))
            final_list.append(entry)

    create_csv(final_list)
    

