# FASTA Sequence Analyzer

A Python command-line tool for analyzing DNA sequences stored in FASTA files. The program calculates basic sequence statistics and performs several common sequence analysis tasks, including transcription, translation, reverse complement generation, six-frame translation, and basic ORF detection.

## Features

* Read DNA sequences from FASTA files
* Validate DNA sequences
* Calculate sequence length
* Calculate GC content
* Count nucleotide composition
* Generate reverse complements
* Transcribe DNA to RNA
* Translate RNA into protein sequences
* Perform translation across the three forward reading frames
* Perform translation across the three reverse-complement reading frames
* Find basic open reading frames (ORFs)
* Translate identified ORFs into proteins
* Calculate amino acid composition
* Identify the longest ORF and corresponding protein
* Analyze multiple sequences in a single FASTA file
* Generate a summary table using Pandas
* Run from the command line using `argparse`
* Test the main functions with Pytest

## Biological Concepts

This project applies basic molecular biology concepts to DNA sequence analysis.

### GC Content

GC content is the percentage of bases in a DNA sequence that are either guanine (G) or cytosine (C).

```text
(G + C) / total bases × 100
```

GC content provides a basic way to describe and compare DNA sequences. Different genomic regions and organisms can have different GC content distributions.

### Nucleotide Composition

Nucleotide composition counts the number of A, T, G, and C bases in a sequence. This provides a simple description of the sequence and can be used alongside GC content when comparing sequences.

### Reverse Complement

The reverse complement represents the complementary DNA strand in the opposite direction. A pairs with T, while G pairs with C.

Reverse complements are important because DNA is double-stranded and biological sequences can be analyzed from either strand.

### Transcription

Transcription converts a DNA sequence into an RNA sequence. In this project, transcription is represented by replacing thymine (T) with uracil (U).

For example:

```text
DNA:  ATGC
RNA:  AUGC
```

### Translation

Translation converts an RNA sequence into a protein sequence using groups of three nucleotides called codons. Each codon corresponds to an amino acid, while specific codons serve as start or stop signals.

This project uses a codon table to translate RNA sequences and stops translation when a stop codon is encountered.

### Reading Frames

A reading frame determines where a nucleotide sequence is divided into groups of three for translation.

A DNA sequence has three possible forward reading frames depending on where the grouping begins. The reverse-complement strand has three additional reading frames, giving six possible reading frames in total.

The project translates all six frames to demonstrate how changing the reading frame can produce different protein sequences.

### Open Reading Frames

An open reading frame (ORF) is a sequence that begins with a start codon and continues until a stop codon is reached within the same reading frame.

This project identifies basic ORFs beginning with the RNA start codon `AUG` and ending at one of the stop codons:

```text
UAA
UAG
UGA
```

Identified ORFs are then translated into protein sequences, and the longest detected ORF is reported.

ORF detection is a basic approach to identifying regions that could potentially encode proteins. In real genomic analysis, additional biological evidence and more sophisticated methods are generally needed to determine whether an ORF represents a functional gene.

## Technologies

* Python
* Pandas
* Pytest
* argparse

## Project Structure

```text
fasta-sequence-analyzer/
├── src/
│   └── analyzer.py
├── tests/
│   └── test_analyzer.py
├── data/
│   └── example.fasta
├── results/
├── requirements.txt
├── LICENSE
└── README.md
```

## Installation

Clone the repository:

```bash
git clone https://github.com/nathanperkins-dev/fasta-sequence-analyzer.git
```

Move into the project directory:

```bash
cd fasta-sequence-analyzer
```

Create a virtual environment:

```bash
python3 -m venv .venv
```

Activate the virtual environment:

```bash
source .venv/bin/activate
```

Install the required packages:

```bash
pip install -r requirements.txt
```

## Usage

Run the analyzer by providing a FASTA file as a command-line argument:

```bash
python src/analyzer.py data/example.fasta
```

The program analyzes each valid DNA sequence and prints the results to the terminal.

For each sequence, the program reports:

* Sequence length
* GC content
* Nucleotide composition
* Reverse complement
* RNA sequence
* Forward reading-frame translations
* Reverse-complement reading-frame translations
* Full translated protein sequence
* Detected ORFs
* ORF protein sequences
* Number of ORFs
* Longest ORF
* Longest protein
* Amino acid composition

A Pandas summary table is also generated for the valid sequences.

## Example

The included `data/example.fasta` file contains three sequences. One contains an invalid nucleotide (`X`), which demonstrates the sequence validation feature.

```text
>sequence_1
ATGCGTACGTAGATGAAATTTTAG
>sequence_2
ATGCGTACGTAGXATGAAATTTTAG
>sequence_3
ATGAAATTTGGGCCCTAG
```

Run the analyzer with:

```bash
python src/analyzer.py data/example.fasta
```

The program identifies `sequence_2` as invalid and skips it during the analysis.

The valid sequences are then included in the Pandas summary:

```text
Sequence Summary:
         name  length     gc  A  T  G  C  number_of_orfs  longest_orf_length  longest_protein_length
0  sequence_1      24  33.33  8  8  6  2               2                  12                       3
1  sequence_3      18  44.44  5  5  5  3               1                  18                       5
```

The program also reports overall statistics such as the longest sequence, average sequence length, average GC content, and total number of detected ORFs.

## Testing

The project includes a Pytest test suite covering the main sequence analysis functions as well as several edge cases.

Run the full test suite with:

```bash
pytest
```

The current test suite contains 18 tests covering:

* Sequence length
* GC content
* Empty-sequence GC handling
* Nucleotide composition
* Reverse complement
* Sequence validation
* Transcription
* Translation
* Translation stopping at stop codons
* Incomplete codons
* Reading-frame translation
* Amino acid composition
* ORF detection
* Sequences without ORFs
* Sequences with a start codon but no stop codon
* Full sequence analysis
* Sequence analysis when no ORFs are present

Current test result:

```text
18 passed
```

## What I Practiced

This project was built to practice applying Python programming to basic bioinformatics problems.

Some of the main programming concepts used include:

* Functions
* Loops
* Conditional statements
* Dictionaries
* Lists and tuples
* String manipulation
* File handling
* Command-line arguments
* Data validation
* Pandas DataFrames
* Unit testing with Pytest
* Git and GitHub

The project also provided practice with breaking a larger problem into smaller functions and then combining those functions into a working command-line program.

On the biology side, the project provided practice working with DNA and RNA sequences, nucleotide composition, GC content, reverse complements, transcription, translation, reading frames, and basic ORF identification.
