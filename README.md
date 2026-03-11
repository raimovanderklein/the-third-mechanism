# The Third Mechanism

Reproducible analysis code for:

**"The Third Mechanism: Evidence from the Genetic Code for a Generative Property of Dissipative Cycling"**

Raimo van der Klein (2026)

Paper: [Zenodo DOI forthcoming]

---

## What this is

The genetic code is organised by two operations — structural stabilisation and functional interaction — in a mandatory cyclic order. This repository contains the code that reproduces every quantitative claim in the paper.

## Scripts

| Script | What it tests | Key result |
|---|---|---|
| `code/genetic_code.py` | Full analysis: hydropathy by codon position, four functional groups, threshold locations, TP53 domain mapping | F = 59.0, p = 1.77 × 10⁻¹⁷ |
| `code/invariance_test.py` | Whether the hold/cross axis is preserved across all 37 known genetic code variants | A-axis: zero violations in 4 billion years |
| `code/mandatory_order.py` | Whether G → T → A → C is the only valid cyclic order out of 24 possible arrangements | Exactly one cycle satisfies all four chemical constraints |
| `code/epigenetic_drift.py` | Whether epigenetic drift correlates with cell division rate, not chronological time | Spearman r = 0.996, p < 0.0001 |

## Requirements

```
Python 3.8+
scipy
matplotlib
numpy
```

Install:
```bash
pip install scipy matplotlib numpy
```

## Running

```bash
cd code
python genetic_code.py       # Full analysis + plots
python invariance_test.py    # Invariance across all life
python mandatory_order.py    # Mandatory cyclic order proof
python epigenetic_drift.py   # Epigenetic drift correlation
```

All scripts print their results to the terminal. `genetic_code.py` also generates plots in `results/`.

## Results

The `results/` directory contains pre-generated outputs:

- `genetic_code_hold_cross.png` — Codon table coloured by hydropathy (Kyte-Doolittle)
- `64_codons_cycle.png` — All 64 codons arranged in cycle order
- `codon_usage_along_gene.png` — Codon regime usage along TP53
- `codon_analysis.csv` — Raw data: every codon with its amino acid, hydropathy, and regime assignment

## Companion papers

- van der Klein, R. (2026). *The Code Truth.* Zenodo. DOI: [10.5281/zenodo.18929161](https://doi.org/10.5281/zenodo.18929161)
- van der Klein, R. (2026). *There Is Only One Way to Grow.* Zenodo. DOI: [10.5281/zenodo.18929428](https://doi.org/10.5281/zenodo.18929428)

## License

MIT
