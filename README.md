# Coupled_Canopy: Solve coupled *C*<sub>i</sub>, *g*<sub>s</sub>, *E*, *A*<sub>n</sub> and *T*<sub>leaf</sub>.

## Overview

Iteratively solve the intercellular CO<sub>2</sub> concentration (*C*<sub>i</sub>), stomatal conductance (*g*<sub>s</sub>), transpiration (*E*), net leaf assimilation (*A*<sub>n</sub>) and leaf temperature (*T*<sub>leaf</sub>).


## Structure

```plaintext
coupled_canopy/
├── __init__.py
├── models/              # Model code
│   ├── __init__.py
│   ├── farquhar.py
│   └── penman_monteith_leaf.py
├── utils/               # Utility functions and constants
│   ├── __init__.py
│   ├── constants.py
│   └── utils.py
examples/                # Example scripts 
├── plot_A_vs_Ca.py
└── ...

## Running examples

From the top directory, run:

```bash
python -m examples.plot_A_vs_Ca

## Contacts

- Martin De Kauwe: mdekauwe at gmail.com
