import sys
sys.path.append("src")

from analyzer import (
    sequence_length,
    gc_content,
    nucleotide_composition,
    reverse_complement,
    validate_sequence,
    transcribe,
    translate,
    translate_frame,
    amino_acid_composition,
    find_orfs,
    analyze_sequence
)

def test_sequence_length():
    assert sequence_length("ATGC") == 4

def test_gc_content():
    assert gc_content("ATGC") == 50.0

def test_gc_content_empty():
    assert gc_content("") == 0.0

def test_nucleotide_composition():
    assert nucleotide_composition("AATGCC") == {
        "A": 2,
        "T": 1,
        "G": 1,
        "C": 2
    }

def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"

def test_validate_sequence():
    assert validate_sequence("ATGC") is True
    assert validate_sequence("ATGX") is False

def test_transcribe():
    assert transcribe("ATGC") == "AUGC"

def test_translate():
    assert translate("AUGCGUUAA") == "MR"

def test_translate_frame():
    assert translate_frame("AUGCGUUAA", 0) == "MR"

def test_translate_stops_at_stop_codon():
    assert translate("AUGUAAAUG") == "M"

def test_amino_acid_composition():
    assert amino_acid_composition("MRTM") == {
        "M": 2,
        "R": 1,
        "T": 1
    }

def test_find_orfs():
    assert find_orfs("AUGCGUUAA") == [
        (0, 9, "AUGCGUUAA")
    ]

def test_analyze_sequence():
    results = analyze_sequence("ATGC")

    assert results["length"] == 4
    assert results["gc"] == 50.0
    assert results["composition"] == {
        "A": 1,
        "T": 1,
        "G": 1,
        "C": 1
    }
    assert results["reverse"] == "GCAT"
    assert results["rna"] == "AUGC"

def test_translate_empty():
    assert translate("") == ""

def test_translate_incomplete_codon():
    assert translate("AUGC") == "M"

def test_find_orfs_no_orf():
    assert find_orfs("CCCGGGCCCAAA") == []

def test_find_orfs_no_stop_codon():
    assert find_orfs("AUGCGUACG") == []

def test_analyze_sequence_no_orf():
    results = analyze_sequence("CCCGGGCCCAAA")

    assert results["orfs"] == []
    assert results["orf_proteins"] == []
    assert results["longest_orf"] is None
    assert results["longest_protein"] is None

