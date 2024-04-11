import csv

def profile_most(text, k, profile):
    max = float('-inf')
    for i in range(len(text)-(k-1)):
        kmer = text[i:i+k]
        score = 1
        for j in range(k):
            cur = kmer[j]
            if cur == "A":
                score *= profile[0][j]
            elif cur == "C":
                score *= profile[1][j]
            elif cur == "G":
                score *= profile[2][j]
            else:
                score *= profile[3][j]
        if score > max:
            max = score
            res = kmer
    return res

def motifmatrix(dna, k):
    matrix = []
    for string in dna:
        row = []
        for i in range(k):
            row.append(string[i])
        matrix.append(row)
    return matrix

def count_from(motifs):
    count = [[0 for j in range(len(motifs[0]))] for i in range(4)]
    for j in range(len(motifs[0])):
        for i in range(len(motifs)):
            cur = motifs[i][j]
            if cur == "A":
                count[0][j] += 1
            elif cur == "C":
                count[1][j] += 1
            elif cur == "G":
                count[2][j] += 1
            else:
                count[3][j] += 1
    return count
        
        

def profile_from(motifs):
    count = count_from(motifs)
    profile = [[0 for j in range(len(count[0]))] for i in range(4)]
    rows = len(motifs)
    for j in range(len(count[0])):
        for i in range(4):
            profile[i][j] = count[i][j] / rows
    return profile

def score(motifs):
    score = 0
    rows = len(motifs)
    for j in range(len(motifs[0])):
        occur = [0 for i in range(4)]
        for i in range(rows):
            curr = motifs[i][j]
            if curr == "A":
                occur[0] += 1
            elif curr == "C":
                occur[1] += 1
            elif curr == "G":
                occur[2] += 1
            else:
                occur[3] += 1
        score += rows - max(occur)
    return score
# Insert your greedy_motif_search function here, along with any subroutines you need
def greedy_motif_search(dna: list[str], k: int, t: int) -> list[str]:
    """Implements the GreedyMotifSearch algorithm."""
    bestmotifs = motifmatrix(dna, k)
    for i in range(len(dna[0])-(k-1)):
        motifs = []
        motifs.append(dna[0][i:i+k])
        for j in range(1, t):
            # print(dna[j])
            profile = profile_from(motifs[0:j])
            motifs.append(profile_most(dna[j], k, profile))
        if score(motifs) <= score(bestmotifs):
            bestmotifs = motifs
    return bestmotifs

def read_fasta(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()

    reads = []
    s = ""
    for read in lines:
        if read[0] == '>':
            continue
        read = read.strip()
        s = s + read
        if len(read) == 1:
            s = s[85:115]
            reads.append(s)
            s = ""
    return reads

def read_full_fasta(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()

    reads = []
    s = ""
    for read in lines:
        if read[0] == '>':
            continue
        read = read.strip()
        s = s + read
        if len(read) == 1:
            reads.append(s)
            s = ""
    return reads

def create_pwm(list_of_kmers, k): 
    pwm = []
    denom = len(list_of_kmers) + 4
    for i in range(k):
        a_count = 1
        c_count = 1
        t_count = 1
        g_count = 1
        for kmer in list_of_kmers:
            match kmer[i]:
                case "A":
                    a_count += 1
                case "C":
                    c_count += 1
                case "G":
                    g_count += 1
                case "T":
                    t_count += 1
        probabilities = {"A": a_count/denom, "C": c_count/denom, "G": g_count/denom, "T": t_count/denom}
        pwm.append(probabilities)
    return pwm


def create_csv(values):
      with open("predictions.csv", 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        for value in values:
            csv_writer.writerow([value])

if __name__ == "__main__":
    k = 13
    reads = read_fasta("boundcentered.fasta")
    test_reads = read_full_fasta("boundcentered.fasta")
    output = greedy_motif_search(reads, k, 1500)
    # print(output)
    pwm = create_pwm(output, k)
    # print(pwm)
    full_reads = read_full_fasta("boundrandomoffset.fasta")

    read_to_use = full_reads

    final = []
    just_nums = []
    for i in range(len(read_to_use)):
        read = read_to_use[i]
        max = 0
        pos = 31
        for j in range(len(read)-k):
            score = 1
            for a in range(k):
                score *= pwm[a][read[j+a]]
            if score > max:
                max = score
                pos = j
        if pos < 31:
            pos = 31
        if pos > 176:
            pos = 176
        entry = "seq" + str(i+1) + "\t" + str(pos)
        just_nums.append(pos)
        # print(entry)
        final.append(entry)

    print(sum(just_nums)/len(just_nums))
    create_csv(final)


    
    