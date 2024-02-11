# Metagenomics Project Code 
# Author: Devon Shao  
# Date: 2/10/2024

#See Read.me for project description
def map_reads_to_genomes(genomes, reads, k):
    reads_dict = {}
    num_of_genomes = len(genomes)
    num_of_reads = len(reads)

    # Iterate through reads, constructing a dict of dicts holding the first third and second third of the read sequence, the genome that the read matched with most closely and the smith_waterman score between the read and the section of the genome that the read matched most closely with
    for i in range(num_of_reads):  
        print(f'reads_dict construction: {i} / {num_of_reads}')
        read_sequence = reads[i]  # Pull out the sequence of the read
        first_third = read_sequence[0:k]
        second_third = read_sequence[k:2*k]
        reads_dict[read_sequence] = {'first_third':first_third, 'second_third':second_third, 'max_score':-999, 'genome_match':0}  

    # Iterate through each genome, creating a hashmap for each genome where each 16-mer section of the genome is a key, and the position of its first nucleotide is the value. 
        # Then iterate through each read. Index the genome hashmap with the first 16 characters of each read to find the matching position in the genome, and use smith_waterman to calculate how well the overall read matches to that genome section. Do the samew with the second 16 characters of each read. 
    for j in range(num_of_genomes):
        print(f"genome_{j} / {num_of_genomes}")
        genome_sequence = genomes[j]
        genome_hashmap = make_kmer_hashmap(genome_sequence, 16)
        
        # Iterate through each read 
        for k in range(num_of_reads):
            read_sequence = reads[k]
            first_third = reads_dict[read_sequence]['first_third']
            second_third = reads_dict[read_sequence]['second_third']

            if first_third in genome_hashmap:
                positions = genome_hashmap[first_third]

                for position in positions:
                    smith_waterman_score = smith_waterman(read_sequence, genome_sequence[position:position + len(read_sequence)])

                    if smith_waterman_score > reads_dict[read_sequence]['max_score']:
                        reads_dict[read_sequence]['max_score'] = smith_waterman_score
                        reads_dict[read_sequence]['genome_match'] = j
            
            elif second_third in genome_hashmap: 
                positions = genome_hashmap[second_third]

                for position in positions:
                    smith_waterman_score = smith_waterman(read_sequence, genome_sequence[position-k:position-k+len(read_sequence)])

                    if smith_waterman_score > reads_dict[read_sequence]['max_score']:
                        reads_dict[read_sequence]['max_score'] = smith_waterman_score
                        reads_dict[read_sequence]['genome_match'] = j
    return reads_dict


# Function to construct a hashmap from a genome where the keys are all the k-mers of the genome and values are the positions in the genome where the k-mers begin.
def make_kmer_hashmap(genome, k):
    kmer_hashmap = {}
    for kmer_start_pos in range(len(genome) - k + 1):
        kmer_sequence = genome[kmer_start_pos:kmer_start_pos+k]
        if kmer_sequence not in kmer_hashmap:
            kmer_hashmap[kmer_sequence] = [kmer_start_pos]
        else: 
            kmer_hashmap[kmer_sequence].append(kmer_start_pos)
    return kmer_hashmap

# This function takes two strings and returns a score that quantifies how good of a match those strings make.
def smith_waterman(v:str, w:str):
    s = [[0] * (len(w)+1) for _ in range(len(v)+1)]

    for i in range(1, len(v)+1):
        s[i][0] = 0
    for j in range(1, len(w)+1):
        s[0][j] = 0
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            match = 0
            if v[i-1] == w[j-1]:
                match = 1
        
            s[i][j] = max(s[i-1][j]-1, s[i][j-1]-1, s[i-1][j-1] + match)

    return s[len(v)][len(w)]


# Example use of the map_reads_to_genomes function, where you provide a list of genomes and a list of reads and k = 16

#predictions = map_reads_to_genomes(genomes, reads, 16)
