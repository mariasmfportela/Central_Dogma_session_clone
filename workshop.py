# My solution to the tasks
# 9th December 2020

# Part 2

# Function to convert a FASTA file into a DNA sequence string object
def fasta_to_string(filename):
    # filename is a string specifying the name of a .fasta file

    # Open file - Part 1
    fasta_file = open(filename)

    # Create empty string
    string = ''

    # Read each line of the FASTA file
    for line in fasta_file:

        # Ignore lines starting in ">"
        if line[0] == '>':
            pass

        # Add characters to the string
        else:
            line = line[:-1]
            string += line

    # Return FASTA sequence as a string
    return string


# Part 3

# Codon to amino acid reference table
codon_dict = {  # from https://www.geeksforgeeks.org/dna-protein-python-3/
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}


# Function to translate a given DNA sequence into a peptide sequence
def mrna_to_protein(mrna_string, codon_start):
    # mrna_string is a string containing a DNA sequence
    # codon_start specifies the first base in the intended reading frame (0-2)

    # Initialize empty protein sequence string
    protein_sequence = ''

    # Iterate through given sequence and translate into protein
    while codon_start < len(mrna_string) - 2:
        # Define the current codon
        codon = mrna_string[codon_start:codon_start + 3]

        # Translate into amino acid
        amino_acid = codon_dict[codon]

        # Add amino acid to protein sequence
        protein_sequence += amino_acid

        # A more accurate way to perform translation, taking start and stop codons into account
        # This exercise files uses an arbitrary motif - so this version won't find it!
        # if amino_acid == 'M':
        #     # This is the start codon so always add
        #     protein_sequence += amino_acid
        #
        # elif amino_acid == '*':
        #     # This is the stop codon so terminate translation
        #     break
        #
        # elif len(protein_sequence) > 0:
        #     # Translation already started so add this amino acid
        #     protein_sequence += amino_acid
        #
        # else:
        #     # Translation hasn't started so skip this amino acid
        #     pass

        # Increment codon start position
        codon_start += 3

    # Return protein sequence as a sequence
    return protein_sequence


# Part 4

# Complement nucleotide base reference table
complement_dict = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'
}


# Function to generate the (reverse) complementary strand from a given DNA sequence
def dna_to_reverse_complement(dna):
    # dna is a string containing a DNA sequence

    # Initialize an empty reverse complement string
    reverse_complement = ''

    # Initialize current base at the end of the sequence
    current_base_index = len(dna) - 1

    # Iterate from the end of the input dna sequence and find complement for each base
    while current_base_index >= 0:
        # Define current base
        base = dna[current_base_index]

        # Determine the complement of the current base
        base_complement = complement_dict[base]

        # Add the base complement to the reverse complement sequence
        reverse_complement += base_complement

        # Decrement curent base index
        current_base_index -= 1

    # Return the reverse complement as a string
    return reverse_complement


# Function to translate all reading frames (RFs) from a given DNA sequence
def get_all_rf(dna):
    # dna is a string containing a DNA sequence

    # Determine the reverse complement of the sequence
    dna_reverse = dna_to_reverse_complement(dna)

    # Create a dictionary of reading frames
    rf_list = {
        'Original 1': mrna_to_protein(dna, 0),
        'Original 2': mrna_to_protein(dna, 1),
        'Original 3': mrna_to_protein(dna, 2),
        'Reverse 1': mrna_to_protein(dna_reverse, 0),
        'Reverse 2': mrna_to_protein(dna_reverse, 1),
        'Reverse 3': mrna_to_protein(dna_reverse, 2),
    }

    # Return reading frames
    return rf_list


# Function calls for parts 1 to 4
spike_protein_mrna = fasta_to_string('sars-cov2_spike_protein_mRNA.fasta')
spike_protein_rf = get_all_rf(spike_protein_mrna)
print(spike_protein_rf)


# Part 5

# Function to search a given motif in a protein sequence
def search_motif(motif, protein):
    # Both motif  and protein are strings representing amino acid sequences

    # Define search initial search index
    index = 0

    # Initialize empty hits object
    hits = []

    # Iterate through protein sequence to find motif
    while index + len(motif) <= len(protein):
        query = protein[index:(index+len(motif))]
        if query == motif:
            hits.append({'start': index + 1, 'end': index + len(motif)})

        index += 1

    return hits


# Part 6

# Function to search for target motif in all reading frames (rf) and print result to console
def search_motif_all_rf(motif, rf_dict):
    # motif is a string representing a peptide sequence
    # rf_dict is a dictionary of peptide sequences for each rf

    # Iterate over all rf
    for rf in rf_dict:
        # Search motif in current rf
        results = search_motif(motif, rf_dict[rf])

        # Print each result
        for i in range(len(results)):
            result = results[i]
            print("Motif start: " + str(result['start']) + "\nMotif end: "
                  + str(result['end']) + "\nReference frame: " + rf)


target_motif = fasta_to_string('sars-cov2_arbitrary_motif.fasta')
print(target_motif)
print("Find arbitrary motif in spike protein")
search_motif_all_rf(target_motif, spike_protein_rf)

# Part 7
a2d_mrna = fasta_to_string('a2d_mRNA.fasta')
a2d_all_rf = get_all_rf(a2d_mrna)

vwa = fasta_to_string('VWA_domain.fasta')
vwa_like = fasta_to_string('VWA-like_domain.fasta')
cache = fasta_to_string('cache_domain.fasta')

print("\nFind motifs in alpha-2-delta-1 protein")
print("\nVWA motif")
search_motif_all_rf(vwa, a2d_all_rf)

print("\nVWA-like motif")
search_motif_all_rf(vwa_like, a2d_all_rf)

print("\nCache motif")
search_motif_all_rf(cache, a2d_all_rf)
