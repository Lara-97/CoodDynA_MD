**Version:** 1.0.0

# Aggregate Analyzer

**Aggregate Analyzer** is a Python-based tool for analyzing molecular simulation data. It detects aggregates from trajectory files, identifies the number of residues involved, and generates plots.

---

## Project Structure
your-project/
├── analyse_aggregates.py	# Main entry point for running analysis
├── config.yaml			# Input configuration
├── requirements.txt 		# Python dependencies
├── outputs/ 			# Output directory (auto-created)
├── src/ 			# Core source code (functions and classes)
│ ├── main.py
│ ├── analysis.py
│ └── plotting.py
├── examples/ 			# Usage examples
│ └── example1/
│ ├── config.yaml
│ └── run_analysis.py

---

## Getting Started

### 1. Clone the Repository
git clone https://github.com/yourusername/aggregate-analyzer.git
cd aggregate-analyzer

### 2. Install Dependencies
pip install -r requirements.txt

### 3. Run the Analysis
python aggregate_analysis_v1_0.py

Or use one of the examples:

cd examples/example1
python run_analysis.py
