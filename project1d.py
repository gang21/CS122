import csv
import math
from collections import defaultdict

# read reference genome
# read pair reads
def read_paired_reads():
    with open("project1d_reads.fasta", 'r') as file:
        reads = file.readlines()
    return reads

def read_ref_genome(filename):
    with open(filename, 'r') as file:
            ref_genome = ''.join([line.strip() for line in file.readlines()[1:]])
    return ref_genome

def smith_waterman(v, w):
    n = len(v)
    m = len(w)
    s = [[0 for i in range(m+1)] for j in range(n+1)]
    Backtrack = [[0 for i in range(m+1)] for j in range(n+1)]
    for i in range(n+1):
        s[i][0] = 0
    for j in range(m+1):
        s[0][j] = 0
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = 0
            if v[i-1] == w[j-1]:
                match = 1
            val = max(s[i-1][j], s[i][j-1], s[i-1][j-1]+match)
            s[i][j] = val
    return s[i][j]

def create_hash_table(reference_genome, kmer_length):
    n = len(reference_genome) - kmer_length + 1
    hash_table = defaultdict(list)
    for i in range(n):
        kmer = reference_genome[i:i+kmer_length]
        hash_table[kmer].append(i)
    return hash_table

def create_csv(values):
      with open("predictions.csv", 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        for value in values:
            csv_writer.writerow([value])

def main():
    paired_reads = read_paired_reads()
    print(paired_reads)
    n = len(paired_reads)
    num_genomes = 5000 # 10 for sample 5000 for 1d
    maxes_per_read = {}
    for i in range(num_genomes):
        # filename = "/Users/gabriellaang/Desktop/CS CM122/Project1/project1c_sample/project1c_sample_genome_" + str(i) + ".fasta"
        filename = "project1d_genome_" + str(i) + ".fasta"
        ref_genome = read_ref_genome(filename)
        hash_table = create_hash_table(ref_genome, 16)
        for j in range(1, n, 2): # check every read with this paired
            read = paired_reads[j]
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
            # check smith-waterman for all positions
            maxes = []
            for position in potential_kmers:
                genome = ref_genome[position:position+len_read]
                sw = smith_waterman(genome, read)
                maxes.append({sw:position})
            if len(maxes) != 0:
                try:
                    maximum = max(maxes)
                except:
                    maximum = maxes[0]
                print("Genome: ", i, "  |  Read: ", math.floor(j/2), " -- ", list(maximum.values())[0])
                key = math.floor(j/2)
                max_val = list(maximum.keys())[0]
                max_pos = list(maximum.values())[0]
                if key in maxes_per_read.keys():
                    read_dict = maxes_per_read[key]
                    get_max = read_dict['max']
                    if max_val > get_max:
                        value = {"reference_genome": i, "max": max_val, "position": max_pos}
                        maxes_per_read[key] = value
                else:
                    # print("we r in the else")
                    value = {"reference_genome": i, "max": max_val, "position": max_pos}
                    maxes_per_read[key] = value

    outputs = []
    for i in range(math.floor(n/2)):
        if i in maxes_per_read.keys():
            print("READ ", str(i))
            output = ">read_" + str(i) + " Genome_Number_" + str(maxes_per_read[i]['reference_genome'])
            outputs.append(output)
        else:
            output = ">read_" + str(i) + " Genome_Number_0"
            outputs.append(output)
    print(outputs)
    create_csv(outputs)



if __name__ == "__main__":
    main()
