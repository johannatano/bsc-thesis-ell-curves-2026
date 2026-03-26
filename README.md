# bsc-thesis-ell-curves-2026

Tools for enumerating elliptic curves over finite fields, organizing them into
isogeny classes, and experimenting with elliptic-point counts and Hecke-style
computations.

## Requirements

- SageMath
- Python packages used by the scripts, in particular `sympy`, `numpy`,
	`requests`, and `tqdm` when progress bars are enabled

Run the command-line tools with Sage's Python:

- `sage -python classifier.py ...`
- `sage -python hecke.py ...`

## Classifier

The main entry point for curve enumeration is [classifier.py](classifier.py).

It supports three practical modes:

- direct enumeration of `j`-invariants and their twists,
- HCP/CM enumeration via Hilbert class polynomials,
- class-number mode for more aggregated counting.

### Basic usage

Run over a specific field $F_{p^n}$ and a specific level $\ell$:

- `sage -python classifier.py -p 7 -n 2 -l 3`

Use HCP/CM enumeration instead of direct geometric enumeration:

- `sage -python classifier.py -p 7 -n 2 -l 3 --use-hcp`

Use class-number mode:

- `sage -python classifier.py -p 7 -n 2 -l 3 --use-cn`

Scan several levels automatically:

- `sage -python classifier.py -p 7 -n 2`

Scan default primes and default extension degrees:

- `sage -python classifier.py`

### Important options

- `-p`, `--p`: field characteristic
- `-n`, `--n`: extension degree
- `-l`, `--l`: level $\ell$
- `--use-hcp`: use Hilbert class polynomial enumeration
- `--use-cn`: use class-number counting instead of explicit curves where possible
- `--rank-method {auto,div_poly,mod_poly,invariants}`: choose how above-floor
	rank tests are performed
- `--true-height`: use exact BFS height in the volcano instead of the default
	floor test

### What it prints

For each selected field and level, the classifier prints the total count of
elliptic-point pairs `(E, P)` at level $\ell$.

## Hecke script

Entry point is [hecke.py](hecke.py).

### Basic usage

Compute at one prime and one level:

- `sage -python hecke.py -p 11 -l 3 -k 2`

Use HCP/CM setup:

- `sage -python hecke.py -p 11 -l 3 -k 4 --use-hcp`

Scan several levels automatically:

- `sage -python hecke.py -p 11 -k 2`

### Important options

- `-p`, `--p`: field characteristic
- `-l`, `--l`: level $\ell$
- `-k`, `--k`: weight
- `--use-hcp`: use Hilbert class polynomial setup
- `--rank-method {auto,div_poly,mod_poly,invariants}`
- `--true-height`

## Repository layout

- [classifier.py](classifier.py): command-line script for curve classification
- [hecke.py](hecke.py): command-line script for Hecke-style experiments
- [lib/](lib): custom lib for number field and curve classifications
- [utils/](utils): shared helpers and older utility code