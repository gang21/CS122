from collections import Counter
import random
import copy
import csv

with open('test.fasta', 'r') as f:
    lines = f.readlines()

sequences = []
sequence = ''
for line in lines:
    if line.startswith('>'):  
        if sequence:  
            sequences.append(sequence)
            sequence = ''
    else: 
        sequence += line.strip()

sequences.append(sequence)

dna = []

with open('bound.fasta', 'r') as f:
    dnalines = f.readlines()
    

string = ''
for line in dnalines:
    if line.startswith('>'):  
        if string:  
            dna.append(string[85:115])
            string = ''
    else: 
        string += line.strip()

dna.append(string)

# with open('notbound.fasta', 'r') as f:
#     dnalines = f.readlines()
    
# unbound = []
# string = ''
# for line in dnalines:
#     if line.startswith('>'):  
#         if string:  
#             dna.append(string)
#             string = ''
#     else: 
#         string += line.strip()

# unbound.append(string)

def reverse_complement(text):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_text = text[::-1]
    reverse_complement_text = ''.join([complement[n] for n in reverse_text])
    return reverse_complement_text

def motifmatrix(dna, k):
    matrix = []
    for string in dna:
        r = random.randint(0, len(string) - k)
        row = []
        for i in range(r, r+k):
            row.append(string[i])
        matrix.append(row)
    return matrix

def count_from(motifs):
    nucleotides = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    count = [[1 for j in range(len(motifs[0]))] for i in range(4)]
    for j in range(len(motifs[0])):
        for i in range(len(motifs)):
            cur = motifs[i][j]
            count[nucleotides[cur]][j] += 1
    return count
        
def profile_from(motifs):
    count = count_from(motifs)
    rows = len(motifs)
    profile = [[count[i][j] / (rows*2) for j in range(len(count[0]))] for i in range(4)]
    return profile

def prob_most_kmer(text, k, profile):
    nucleotides = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    scores = []
    for i in range(len(text)-(k-1)):
        kmer = text[i:i+k]
        score = 1
        for j in range(k):
            cur = kmer[j]
            score *= profile[nucleotides[cur]][j]
        scores.append(score)

    # Compute the sum of all scores to normalize them to probabilities
    total_score = sum(scores)
    probabilities = [score / total_score for score in scores]

    # Use Roulette Wheel Selection to choose a kmer according to its score probability
    random_num = random.uniform(0, 1)
    probability_sum = 0
    for i in range(len(probabilities)):
        probability_sum += probabilities[i]
        if probability_sum >= random_num:
            result = text[i:i+k]
            break

    return result

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
    reverse = reverse_complement(text)
    for i in range(len(reverse)-(k-1)):
        kmer = reverse[i:i+k]
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
    return max

def hamming_distance(p,q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count

def score(motifs):
    score = 0
    rows = len(motifs)
    for j in range(len(motifs[0])):
        occur = Counter(motifs[i][j] for i in range(rows))
        score += rows - occur.most_common(1)[0][1]
    return score

def greedymotifsearch(dna, k, t):
    motifs = motifmatrix(dna, k)
    bestmotifs = motifs
    iter = 0
    for j in range(0, t):
        i = random.randint(0,len(dna)-1)
        motifs.remove(motifs[i])
        profile = profile_from(motifs)
        motifs.insert(i, prob_most_kmer(dna[i], k, profile))
        if score(motifs) < score(bestmotifs):
            bestmotifs = motifs
            iter = 0
        else:
            iter += 1
        if iter == 10:
            return profile_from(bestmotifs)
    return profile_from(bestmotifs)

def prediction(sequences, pwm):
    scores = {}
    for i in range(len(sequences)):
        scores[i] = profile_most(sequences[i], 23, pwm)
    sorted_scores = dict(sorted(scores.items(), key=lambda x: x[1], reverse=True))
    top_6000_scores = list(sorted_scores.keys())[:6000]
    return top_6000_scores


def create_csv(values):
      with open("predictions.csv", 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        for value in values:
            csv_writer.writerow([value])


num = 20000
print(num)
pwm = greedymotifsearch(dna, 20, num)
top_scores = prediction(sequences, pwm)

final = []
for score in top_scores:
    entry = "seq" + str(score + 1)
    final.append(entry)

create_csv(final)


