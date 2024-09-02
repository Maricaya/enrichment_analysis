#!/bin/env python
import json
import gseapy as gp
import os
import shutil
import argparse

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Prepare databases for GSEApy.")
    parser.add_argument("--input", required=True, help="Path to the input database file (GMT or JSON).")
    parser.add_argument("--output", required=True, help="Path to the output GMT file.")
    parser.add_argument("--database", required=True, help="Name of the database.")
    args = parser.parse_args()

    db_path = args.input
    results_path = args.output
    db = args.database

    # if GMT, just copy
    if db_path.lower().endswith('.gmt'):
        shutil.copy(db_path, results_path)
    elif db_path.lower().endswith('.json'):
        # JSON load and save as GMT
        with open(db_path, 'r') as f:
            data = json.load(f)

        with open(results_path, 'w') as f:
            for key, values in data.items():
                f.write(f"{key}\t\t" + "\t".join(values) + "\n")
    else:
        print("Error: Please provide a GMT (*.gmt) or JSON (*.json) database file.")

if __name__ == "__main__":
    main()