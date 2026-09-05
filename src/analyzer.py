import argparse


def read_fasta(filename):
    sequences = {}
    current_name = None

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith(">"):
                current_name = line[1:]
                sequences[current_name] = ""
            else:
                sequences[current_name] += line

    return sequences


def validate_sequence(sequence):
    for base in sequence:
        if base not in ["A", "T", "G", "C"]:
            return False

    return True


def sequence_length(sequence):
    return len(sequence)


def gc_content(sequence):
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


def nucleotide_composition(sequence):
    composition = {}
    composition["A"] = sequence.count("A")
    composition["T"] = sequence.count("T")
    composition["G"] = sequence.count("G")
    composition["C"] = sequence.count("C")

    return composition


def reverse_complement(sequence):
    complements = {"A": "T", "T": "A", "G": "C", "C": "G"}

    result = ""

    for base in sequence:
        result += complements[base]

    result = result[::-1]

    return result


def transcribe(sequence):
    rna = sequence.replace("T", "U")
    return rna


codon_table = {
    "UUU": "F", "UUC": "F",
    "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y",
    "UGU": "C", "UGC": "C",
    "UGG": "W",

    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "AUU": "I", "AUC": "I", "AUA": "I",
    "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",

    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


def translate(rna):
    protein = ""

    for i in range(0, len(rna), 3):
        codon = rna[i:i + 3]

        if codon in ["UAA", "UAG", "UGA"]:
            break

        amino_acid = codon_table[codon]
        protein += amino_acid

    return protein


def amino_acid_composition(protein):
    composition = {}

    for amino_acid in protein:
        if amino_acid in composition:
            composition[amino_acid] += 1
        else:
            composition[amino_acid] = 1
    return composition


def find_orfs(rna):
    orfs = []

    for start in range(len(rna)):
        if rna[start:start + 3] == "AUG":

            for i in range(start, len(rna), 3):
                codon = rna[i:i + 3]

                if codon in ["UAA", "UAG", "UGA"]:
                    orfs.append(rna[start:i + 3])
                    break

    return orfs



def analyze_sequence(sequence):
    length = sequence_length(sequence)
    gc = gc_content(sequence)
    composition = nucleotide_composition(sequence)
    reverse = reverse_complement(sequence)
    rna = transcribe(sequence)
    protein = translate(rna)
    aa_composition = amino_acid_composition(protein)
    orfs = find_orfs(rna)

    orf_proteins = []

    for orf in orfs:
        orf_protein = translate(orf)
        orf_proteins.append(orf_protein)

    if orfs:
        longest_orf = max(orfs, key=len)
    else:
        longest_orf = None

    if longest_orf:
        longest_orf_index = orfs.index(longest_orf)
        longest_protein = orf_proteins[longest_orf_index]
    else:
        longest_protein = None

    return {
        "length": length,
        "gc": gc,
        "composition": composition,
        "reverse": reverse,
        "rna": rna,
        "protein": protein,
        "aa_composition": aa_composition,
        "orfs": orfs,
        "orf_proteins": orf_proteins,
        "longest_orf": longest_orf,
        "longest_protein": longest_protein
    }




parser = argparse.ArgumentParser(
    description="Analyze DNA and RNA FASTA sequences."
)

parser.add_argument("filename", help="Path to FASTA file")

args = parser.parse_args()

sequences = read_fasta(args.filename)

number_of_sequences = len(sequences)
valid_sequences = 0

total_length = 0
total_gc = 0
total_orfs = 0
longest_sequence = None

for name, seq in sequences.items():
    if not validate_sequence(seq):
        print(f"Invalid DNA sequence: {name}")
        continue

    valid_sequences += 1

    results = analyze_sequence(seq)

    length = results["length"]
    total_length += length

    if longest_sequence is None or length > longest_sequence[1]:
        longest_sequence = (name, length)

    gc = results["gc"]
    total_gc += gc

    composition = results["composition"]
    reverse = results["reverse"]
    rna = results["rna"]
    protein = results["protein"]
    aa_composition = results["aa_composition"]
    orfs = results["orfs"]

    number_of_orfs = len(orfs)
    total_orfs += number_of_orfs

    longest_orf = results["longest_orf"]

    orf_proteins = results["orf_proteins"]

    longest_protein = results["longest_protein"]

    print(f"{name}: {length} bp | GC content: {gc:.2f}%")
    print(f"Nucleotide composition: A={composition['A']} | T={composition['T']} | G={composition['G']} | C={composition['C']}")
    print(f"Reverse complement: {reverse}")
    print(f"RNA sequence: {rna}")
    print(f"Protein sequence: {protein}")
    print(f"ORFs: {orfs}")
    print(f"ORF proteins: {orf_proteins}")
    print(f"Number of ORFs: {number_of_orfs}")

    if longest_orf:
        print(f"Longest ORF: {longest_orf}")
        print(f"Longest ORF length: {len(longest_orf)} nt")
        print(f"Longest protein: {longest_protein}")
        print(f"Longest protein length: {len(longest_protein)} aa")

    for index, orf in enumerate(orfs):
        print(f"ORF {index +1}:")
        print(f" RNA: {orf}")
        print(f" Length: {len(orf)} nt")
        print(f" Protein: {orf_proteins[index]}")
        print(f" Protein length: {len(orf_proteins[index])} aa")

    print(f"Amino acid composition: {aa_composition}")

if valid_sequences > 0:
    average_length = total_length / valid_sequences
    average_gc = total_gc / valid_sequences

    print(f"Longest sequence: {longest_sequence[0]} ({longest_sequence[1]} bp)")
    print(f"Average sequence length: {average_length:.2f} bp")
    print(f"Average GC content: {average_gc:.2f}%")
    print(f"Total ORFs found: {total_orfs}")
else:
    print("No valid DNA sequences found.")

