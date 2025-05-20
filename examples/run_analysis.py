import importlib.util
import os
import sys

# Get the absolute path to the analyse_aggregates.py script (one level up)
script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../aggregate_analysis_v1_0.py"))

# Load analyse_aggregates.py as a module
spec = importlib.util.spec_from_file_location("aggregate_analysis_v1_0", script_path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

# Now call the main function from it
if __name__ == "__main__":
    module.main("config.yaml")
