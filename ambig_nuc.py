# Script created to loop through ambiguous nucleotide combinations
# Created by: Akilia Mathie
# Last Edited: 09/06/24

# Import packages
import argparse
from itertools import product

# Dictionary of Ambiguous IUPAC Nucleotide Code Base Combinations
iupac_dict = {'R':['A','G'],
        'Y':['C','T'],
        'S':['G','C'],
        'W':['A','T'],
        'K':['G','T'],
        'M':['A','C'],
        'B':['C','G','T'],
        'D':['A','G','T'],
        'H':['A','C','T'],
        'V':['A','C','G'],
        'N':['A','C','G','T']
        }

# Function to read the FASTA file and return sequences
def read_fasta(file_path):
    with open(file_path, 'r') as f:
        header = None
        sequence = ''
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:  # If we already have a sequence, yield it
                    yield header, sequence
                header = line
                sequence = ''
            else:
                sequence += line
        yield header, sequence  # Yield the last sequence

# Function to generate all combinations for ambiguous codes
def replace_ambiguous(seq, iupac_dict):
    # Find positions of ambiguous nucleotides
    positions = [i for i, nuc in enumerate(seq) if nuc in iupac_dict]
    if not positions:
        return [seq]  # If no ambiguous nucleotides, return the sequence as is
    
    # Generate all possible nucleotide combinations
    base_combinations = [iupac_dict[seq[i]] for i in positions]
    for replacements in product(*base_combinations):
        # Create a new sequence by replacing ambiguous positions
        new_seq = list(seq)
        for pos, repl in zip(positions, replacements):
            new_seq[pos] = repl
        yield ''.join(new_seq)

# Read args
parser = argparse.ArgumentParser()
parser.add_argument('--file',
                    required=True,
                    help="Fasta file with sequence")
args = parser.parse_args()

try:
    args = parser.parse_args()
except AttributeError:
    parser.print_help()
    exit()
    
# Define path of uploaded file
file = args.file

# Output file to write new sequences
output_file = file.split('.')[0] + "_expanded.fasta"

# Process the FASTA file and generate new sequences
with open(output_file, 'w') as out_f:
    for header, sequence in read_fasta(file):
        out_f.write(f"{header}\n")
        for new_seq in replace_ambiguous(sequence, iupac_dict):
            out_f.write(f"{new_seq}\n") 
