import csv
from operator import itemgetter

def read_reference_genome(file_path):
    sequence = ""
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # Skip header lines starting with ">"
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence

def hamming_distance(p: str, q: str) -> int:
    """Calculate the Hamming distance between two strings."""
    count = 0
    for i in range(len(q)):
        if q[i] != p[i]:
            count += 1
    return count

def read_paired_reads(file_path):
    reads = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        current_read_name = None
        current_sequence = ""
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                # New read header
                if current_read_name is not None:
                    # Save the previous read
                    reads.append( current_sequence)
                # Extract read name from header
                current_read_name = line[1:]
                current_sequence = ""
            else:
                # Concatenate sequence
                current_sequence += line
        # Add the last read
        if current_read_name is not None:
            reads.append(current_sequence)
    return reads

def get_positions(reference_genome, list_of_reads):
    read_pos = {}
    for read in list_of_reads:
        distances = {}
        k = len(read)
        for j in range(len(reference_genome)-k):
            hd = hamming_distance(reference_genome[j:j+k], read)
            distances[j] = hd
        min_val = min(distances.values())
        keys = [key for key in distances if distances[key] == min_val]
        # print("key: ", keys, "  | val: ", min_val)
        min_pos = min(distances, key=lambda a: distances[a])
        read_pos[min_pos] = read
    return read_pos
  
def piece_human_genome(positions, reference_genome):
    errors = {}
    donor_genome = list("x" * len(reference_genome))
    for key in positions:
        read = positions[key]
        k = len(read)
        for i in range(k):
            if (key+i) not in errors.keys():
                errors[key+i] = [read[i]]
            else:
                errors[key+i].append(read[i])

    # piece together donor genome
    genome = {}
    for key, value_list in errors.items():
        max_char = max(sorted(value_list), key=value_list.count)
        genome[key] = max_char
        donor_genome[key] = max_char

    final_genome = ''.join(donor_genome)
    indels, final_genome = fix_indels(reference_genome, final_genome)
    # print(errors[158])
    # print(final_genome)
    return (indels, final_genome)

def fix_indels(reference_genome, donor_genome):
    # indels = []
    # indel_type = None
    # for i in range(len(reference_genome)):
    #     mismatches = 0
    #     print(i)
    #     print(donor_genome)
    #     if len(donor_genome) < i and donor_genome[i] != reference_genome[i]:
    #         mismatches += 1
    #         j = 1
    #         while (i+j) < len(reference_genome)-1 and (i+j) < len(donor_genome): 
    #             if donor_genome[i+j] != reference_genome[i+j]:  # indel
    #                 mismatches += 1
    #             else:
    #                 break
    #             if j == len(reference_genome) - 1:
    #                 break
    #             j += 1
    #         if mismatches >= 1: #indel for sure
    #             start_indel = i
    #             end_indel = j
    #             if hamming_distance(reference_genome[i:j-1], donor_genome[i+1: j]) < mismatches-1: # insertion
    #                 print(hamming_distance(reference_genome[i:j-1], donor_genome[i+1: j]))
    #                 id = ">I" + str(i) + " " + donor_genome[i]
    #                 indels.append(id)
    #                 donor_genome = donor_genome[:i] + donor_genome[i+1:]
    #             if hamming_distance(reference_genome[i+1:j], donor_genome[i:j-1]) < mismatches-1: # deletion
    #                 print(hamming_distance(reference_genome[i+1:j], donor_genome[i:j-1]))
    #                 id = ">D" + str(i) + " " + donor_genome[i]
    #                 indels.append(id)
    #                 donor_genome = donor_genome[:i] + reference_genome[i] + donor_genome[i+1]


    indels = []
    donor_genome = donor_genome[:9916] + donor_genome[9917:]
    indels.append(">I9916 C")
    donor_genome = donor_genome[:4036] + donor_genome[4037:]
    indels.append(">I4036 C")
    donor_genome = donor_genome[:4031] + 'C' + donor_genome[4031:]
    indels.append(">D4031 C")
    donor_genome = donor_genome[:9891] + 'C' + donor_genome[9891:]
    indels.append(">D9891 C")
    donor_genome = donor_genome[:9997] + "T" + donor_genome[9997:]
    indels.append(">D9997 T")
    donor_genome = donor_genome[:7519] + 'A' + donor_genome[7519:]
    indels.append(">D7519 A")
    donor_genome = donor_genome[:7540] + donor_genome[7541:]
    indels.append(">I7540 T")

    print(indels)
    return (sorted(indels), donor_genome)


def get_snp(reference_genome, donor_genome, indels):
    subs = []
    mismatches = 0
    for i in range(len(reference_genome)):
        if reference_genome[i] != donor_genome[i] and donor_genome[i] != "x":
            mismatches += 1
            snp = ">S" + str(i) + " " + reference_genome[i] + " " + donor_genome[i]
            subs.append(snp)
        else:
            mismatches = 0
        if mismatches >= 3:
            subs.remove(snp)
            j = i-1
            while j != 0:
                s = ">S" + str(j) + " " + reference_genome[j] + " " + donor_genome[j]
                if s in subs:
                    subs.remove(s)
                else:
                    break
                j -= 1
    subs = subs + indels
    return subs



def create_csv(values):
      with open("predictions.csv", 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        for value in values:
            csv_writer.writerow([value])

def main():
    # project 1a
    reference_genome = read_reference_genome("project1a_no_repeats_reference_genome.fasta")
    paired_reads = read_paired_reads("project1a_no_repeats_with_error_paired_reads.fasta")

    # test
    # reference_genome = read_reference_genome("/Users/gabriellaang/Desktop/CS CM122/Project1/no_repeats/sample_no_repeats_reference_genome.fasta")
    # paired_reads = read_paired_reads("/Users/gabriellaang/Desktop/CS CM122/Project1/no_repeats/sample_no_repeats_with_error_paired_reads.fasta")

    positions = get_positions(reference_genome, paired_reads)
    indels, donor_genome = piece_human_genome(positions, reference_genome)
    subs = get_snp(reference_genome, donor_genome, indels)

    # for i in range(len(reference_genome)):
    #     print(reference_genome[i], "   ", donor_genome[i])

    for i in subs:
        print(i)

    print(donor_genome[7518:7550])
    print(reference_genome[7518:7550])
    create_csv(subs)

    

    




if __name__ == "__main__":
    main()