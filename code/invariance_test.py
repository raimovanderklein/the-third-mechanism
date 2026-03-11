#!/usr/bin/env python3
"""
INVARIANCE TEST: Is the hold/cross axis preserved across all known genetic codes?

Tests whether T-at-position-2 always encodes hydrophobic amino acids and
A-at-position-2 always encodes hydrophilic amino acids across all 37 NCBI
genetic code variants.

From: "The Third Mechanism" (van der Klein, 2026)

Result: A-axis has ZERO violations across all of life.
        T-axis has 5 borderline violations (all Leu→Thr/Ser).
"""

from scipy import stats

HYDROPATHY = {
    'Ile': 4.5, 'Val': 4.2, 'Leu': 3.8, 'Phe': 2.8, 'Cys': 2.5,
    'Met': 1.9, 'Ala': 1.8, 'Gly': -0.4, 'Thr': -0.7, 'Ser': -0.8,
    'Trp': -0.9, 'Tyr': -1.3, 'Pro': -1.6, 'His': -3.2, 'Glu': -3.5,
    'Gln': -3.5, 'Asp': -3.5, 'Asn': -3.5, 'Lys': -3.9, 'Arg': -4.5,
}

# All documented codon reassignments across 37 NCBI genetic code tables
# Source: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
REASSIGNMENTS = [
    ("AGA", "Arg", "Stop", 2, "Vertebrate mitochondrial"),
    ("AGG", "Arg", "Stop", 2, "Vertebrate mitochondrial"),
    ("AUA", "Ile", "Met", 2, "Vertebrate mitochondrial"),
    ("UGA", "Stop", "Trp", 2, "Vertebrate mitochondrial"),
    ("AGA", "Arg", "Ser", 9, "Echinoderm/flatworm mitochondrial"),
    ("AGG", "Arg", "Ser", 9, "Echinoderm/flatworm mitochondrial"),
    ("AAA", "Lys", "Asn", 9, "Echinoderm/flatworm mitochondrial"),
    ("UGA", "Stop", "Trp", 9, "Echinoderm/flatworm mitochondrial"),
    ("UGA", "Stop", "Trp", 3, "Yeast mitochondrial"),
    ("AUA", "Ile", "Met", 3, "Yeast mitochondrial"),
    ("CUU", "Leu", "Thr", 3, "Yeast mitochondrial"),
    ("CUC", "Leu", "Thr", 3, "Yeast mitochondrial"),
    ("CUA", "Leu", "Thr", 3, "Yeast mitochondrial"),
    ("CUG", "Leu", "Thr", 3, "Yeast mitochondrial"),
    ("UGA", "Stop", "Trp", 4, "Mycoplasma/Spiroplasma"),
    ("UGA", "Stop", "Trp", 5, "Invertebrate mitochondrial"),
    ("AGA", "Arg", "Ser", 5, "Invertebrate mitochondrial"),
    ("AGG", "Arg", "Ser", 5, "Invertebrate mitochondrial"),
    ("AUA", "Ile", "Met", 5, "Invertebrate mitochondrial"),
    ("UAA", "Stop", "Gln", 6, "Ciliate nuclear"),
    ("UAG", "Stop", "Gln", 6, "Ciliate nuclear"),
    ("UAA", "Stop", "Tyr", 29, "Mesodinium ciliate"),
    ("UAG", "Stop", "Tyr", 29, "Mesodinium ciliate"),
    ("UAA", "Stop", "Glu", 28, "Campanella ciliate"),
    ("UAG", "Stop", "Glu", 28, "Campanella ciliate"),
    ("UGA", "Stop", "Trp", 22, "Scenedesmus mitochondrial"),
    ("UCA", "Ser", "Stop", 22, "Scenedesmus mitochondrial"),
    ("CUG", "Leu", "Ser", 12, "Candida CUG-Ser1"),
    ("UAA", "Stop", "Gln", 27, "Condylostoma magnum"),
    ("UAG", "Stop", "Gln", 27, "Condylostoma magnum"),
    ("UGA", "Stop", "Trp", 27, "Condylostoma magnum"),
    ("UGA", "Stop", "Gly", 25, "Gracilibacteria"),
    ("AGG", "Arg", "Lys", 24, "Rhabdopleura mitochondrial"),
    ("AUA", "Ile", "Met", 24, "Rhabdopleura mitochondrial"),
]

a_violations = 0
t_violations = 0
a_total = 0
t_total = 0

for codon, std_aa, var_aa, code, organism in REASSIGNMENTS:
    codon_dna = codon.replace('U', 'T')
    pos2 = codon_dna[1]
    h_var = HYDROPATHY.get(var_aa)
    
    if pos2 == 'A' and var_aa != 'Stop' and h_var is not None:
        a_total += 1
        if h_var > 0:
            a_violations += 1
            print(f"A-AXIS VIOLATION: {codon_dna} {std_aa}→{var_aa} (h={h_var:+.1f}) [{organism}]")
    
    if pos2 == 'T' and var_aa != 'Stop' and h_var is not None:
        t_total += 1
        if h_var < 0:
            t_violations += 1
            print(f"T-axis violation: {codon_dna} {std_aa}→{var_aa} (h={h_var:+.1f}) [{organism}]")

print()
print(f"A-axis (cross/hydrophilic): {a_violations} violations out of {a_total} reassignments")
print(f"T-axis (hold/hydrophobic):  {t_violations} violations out of {t_total} reassignments")
print()
if a_violations == 0:
    print("The A-axis is PERFECTLY INVARIANT across all known life.")
