**Version:** 1.0.0

# Coordination Dynamics Analyzer (CoordDynA_MD)



This repository contains a Python-based analysis tool developed to study the formation and structure of molecular aggregates in molecular dynamics (MD) simulations, specifically targeting ion complexes. The code was designed for high-throughput and robust analysis of cation–ligand coordination, aggregate formation, and structural characterization over time.

🧪 This code was developed as a part of scientific publication.
DOI: [10.1016/j.molliq.2025.128711](https://www.sciencedirect.com/science/article/pii/S0167732225018884)

## What the Code Does

Using the NumPy and MDAnalysis Python libraries, this tool performs the following tasks:

- Identification of coordination spheres:
Detects ligands in the first and second coordination spheres of Eu ions based on distance cutoffs from RDF (radial distribution function) data.
- Local environment filtering:
Reduces computational cost by selecting only atoms within a defined cutoff region for each ion.
- Aggregate detection:
Uses a Union-Find algorithm to identify when ligands bridge multiple ions, indicating the formation of oligomers or larger aggregates.
- Parallelized execution:
Employs Python’s multiprocessing module to analyze multiple ions concurrently, significantly reducing computation time.
- Results output:
Stores data about the number and type of ligands in aggregates, as well as coordination residues, in easily readable files.
- Visualization:
Generates plots showing the aggregate sizes and structures using Matplotlib.

---

## Project Structure

```
micelle_analysis/
├── aggregate_analysis_v1_0.py    # Main entry point for running analysis

├── config.yaml                   # Input data

├── requirements.txt              # Python dependencies

├── outputs/                      # Output directory (auto-created)

├── src/                          # Core source code (functions and classes)

│    ├── main.py

│    ├── analysis.py

│    └── plotting.py

├── examples/                    # Usage examples

│    └── example1/

│    ├── config.yaml

│    └── run_analysis.py
```
---

## Getting Started

### 1. Clone the Repository
git clone https://github.com/Lara-97/CoodDynA_MD/

cd micelle_analysis

### 2. Install Dependencies
pip install -r requirements.txt

### 3. Prepare config.yaml with:
- Paths to input trajectory and topology
- Central ion and ligand selection
- Cutoff distances for ion-ligand and ligand-ligand interactions


### 4. Run the Analysis
python aggregate_analysis_v1_0.py


Or use one of the examples:

cd examples/

python run_analysis.py
