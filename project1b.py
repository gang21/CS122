import csv
from collections import defaultdict

def hamming_distance(p: str, q: str) -> int:
    """Calculate the Hamming distance between two strings."""
    count = 0
    l1 = len(q)
    l2 = len(p)
    for i in range(min(l1, l2)):
        if q[i] != p[i]:
            count += 1
    if l1 != l2:
        return count + abs(l1-l2)
    return count

def read_reference_genome(file_path):
    sequence = ""
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # Skip header lines starting with ">"
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence

def read_paired_reads(file_path):
    with open('/Users/gabriellaang/Desktop/CS CM122/Project1/project1b/project1b-u_with_error_paired_reads.fasta', 'r') as file:
        reads = file.readlines()
    return reads

def create_hash_table(reference_genome, kmer_length):
    n = len(reference_genome) - kmer_length + 1
    hash_table = defaultdict(list)
    for i in range(n):
        kmer = reference_genome[i:i+kmer_length]
        hash_table[kmer].append(i)
    return hash_table

def insertion(read, reference_genome, hash, indels):
    first_positions = hash[read[0:16]]
    last_positions = hash[read[32:48]]
    for first in first_positions:
        for last in last_positions:
            if last - (first+16) == 15:
                pos =  find_insertion(read[16:32], reference_genome[first+16:first+32], first+16)
                insertion = ">I" + str(pos) + " " + read[i]
                return pos

def map(read, reference_genome, hash_table):
    l = len(read)
    indices = []
    perfect_matches = []
    for i in range(0, l-16, 16):
        positions = hash_table[read[i:i+16]]
        perfect_matches.append(positions)
        for position in positions:
            if mismatches(read, reference_genome[position - i:position - i + l]) < 3:
                indices.append(position - i)
    return indices

def mismatches(read, ref):
    if len(read) != len(ref):
        return 3
    mismatches = 0
    for i in range(len(read)):
        if read[i] != ref[i]:
            mismatches += 1
    return mismatches

def key(mutation):
    return int(mutation[2:-3])

def ins_key(ins):
    return int(ins[2:-1])

def del_key(d):
    return int(d[2:-1])

def create_csv(values):
      with open("predictions.csv", 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        for value in values:
            csv_writer.writerow([value])

def find_insertion(read, genome, index):
    for i in range(len(read)):
        if read[i] != genome[i]:
            pos = index + i - 1
            label = ">I" + str(pos) + " " + read[i]
            return label

def insertion(read, reference_genome, hash_table):
    first_positions = hash_table[read[0:16]]
    last_positions = hash_table[read[32:48]]
    for first in first_positions:
        for last in last_positions:
            if last - (first+16) == 15:
                return find_insertion(read[16:32], reference_genome[first+16:first+32], first+16)

def delete(read, reference_genome, hash_table):
    first_positions = hash_table[read[0:16]]
    last_positions = hash_table[read[32:48]]
    for first in first_positions:
        for last in last_positions:
            if last - (first+16) == 17:
                return find_deletion(read[16:32], reference_genome[first+16:first+32], first+16)

def find_deletion(read, genome, index):
    for i in range(len(read)):
        if read[i] != genome[i]:
            pos = index + i + 1
            label = ">D" + str(pos) + " " + genome[i]
            return label

def main():
    reference_genome = read_reference_genome("project1b-u_reference_genome.fasta")
    paired_reads = read_paired_reads("project1b-u_with_error_paired_reads.fasta")

    hash_table = create_hash_table(reference_genome, 16)
    mutations = {}
    insertions = {}
    deletions = {}
    for i in range(1, len(paired_reads), 2):
        read = paired_reads[i].strip()
        indices = map(read, reference_genome, hash_table)
        # print(indices)
        if len(indices) > 0:
            for j in range(len(read)):
                if reference_genome[indices[0]+j] != read[j]:
                    index = indices[0]+j
                    original = reference_genome[indices[0]+j]
                    mutation = read[j]
                    label = ">S" + str(index) + " " + original + " " + mutation
                    if label in mutations:
                        mutations[label] += 1
                    else:
                        mutations[label] = 1
        else:
            ins = insertion(read, reference_genome, hash_table)
            d = delete(read, reference_genome, hash_table)
            # determine whether it's an insertion or deletion
            if ins is not None:
                if ins in insertions:
                    insertions[ins] += 1
                else:
                    insertions[ins] = 1
            elif d is not None:
                if d in deletions: 
                    deletions[d] += 1
                else:
                    deletions[d] = 1
    mutations_list = []
    for label in mutations:
        if mutations[label] >= 2:
            mutations_list.append(label)
    insertions_list = []
    for ins in insertions:
        if insertions[ins] >= 2:
            insertions_list.append(ins)
    deletions_list = []
    for d in deletions:
        if deletions[d] >= 2:
            deletions_list.append(d)
    sorted_mutations = sorted(mutations_list, key=key)
    sorted_insertions = sorted(insertions_list, key=ins_key)
    sorted_deletions = sorted(deletions_list, key = del_key)

    print(sorted_deletions)

    all_mutations = sorted_mutations + sorted_insertions + sorted_deletions
    # print(all_mutations)

    create_csv(all_mutations)

if __name__ == "__main__":
    main()
