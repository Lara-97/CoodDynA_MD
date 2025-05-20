import sys
import os
from pathlib import Path
import yaml

# Step 1: Add the project root (this file's location) to sys.path
project_root = Path(__file__).resolve().parent
sys.path.insert(0, str(project_root))

# Step 2: Import your actual main logic from src/
from src.main import main  # assuming your real logic is in src/main.py

# Step 3: Parse the config file path from command line
def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyse_aggregates.py path/to/config.yaml")
        sys.exit(1)

    config_path = sys.argv[1]
    config = load_config(config_path)

    # Step 4: Run the main analysis logic
    main(config)
