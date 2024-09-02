import os
import subprocess
import yaml
import pandas as pd

# Load the configuration file
config_path = os.path.abspath('test/config/example_enrichment_analysis_config.yaml')
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

# Set result path
result_path = os.path.abspath(config['result_path'])

# Create Conda environment (if it does not exist)
def create_conda_env(env_name, env_file):
    try:
        result = subprocess.run(['conda', 'env', 'list'], capture_output=True, text=True)
        if env_name in result.stdout:
            print(f"Conda environment '{env_name}' already exists.")
        else:
            print(f"Creating Conda environment '{env_name}'...")
            subprocess.run(['conda', 'env', 'create', '-f', env_file], check=True)
            print(f"Conda environment '{env_name}' created successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while creating the Conda environment: {e}")
        exit(1)

# Main function to plot enrichment results
def plot_enrichment_result(feature_set, tool, db, conda_env):
    r_script = 'workflow/scripts/enrichment_plot.R'

    # Log file path
    log_file = os.path.join(result_path, 'logs', f'plot_enrichment_result_{tool}_{feature_set}_{db}.log')

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Print paths for debugging
    print(f"Processing feature set: {feature_set}, tool: {tool}, database: {db}")

    # Construct the command to run the R script
    command = [
        'conda', 'run', '--no-capture-output', '--name', conda_env,
        'Rscript', r_script, feature_set, tool, db
    ]

    # Run the R script
    with open(log_file, 'w') as log:
        result = subprocess.run(command, stdout=log, stderr=log, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Plotting enrichment result failed for feature set '{feature_set}', tool '{tool}', and database '{db}'. Check the log file {log_file} for details.")
    else:
        print(f"Plotting completed successfully for feature set '{feature_set}', tool '{tool}', and database '{db}'.")

# Automatically retrieve feature_set, tool, and db
def get_feature_sets_tools_dbs(config):
    # Load annotation file to get feature_set and tool names
    annotation_path = os.path.abspath(config['annotation'])
    annot = pd.read_csv(annotation_path, index_col='name')

    # Extract available feature_set
    feature_sets = annot.index.unique().tolist()

    # Extract available tools (tool) and databases (db)
    tools = config["column_names"].keys()
    databases = config["local_databases"].keys()

    return feature_sets, tools, databases

if __name__ == "__main__":
    # Create Conda environment
    conda_env = 'visualization'
    env_file = os.path.abspath('workflow/envs/visualization.yaml')
    create_conda_env(conda_env, env_file)

    # Automatically get available feature_set, tool, and db
    feature_sets, tools, databases = get_feature_sets_tools_dbs(config)

    print("feature_sets", feature_sets)
    print("tools", tools)
    print("databases", databases)

    # Loop through all combinations of feature_set, tool, and db
    for feature_set in feature_sets:
        for tool in tools:
            for db in databases:
                try:
                    plot_enrichment_result(feature_set, tool, db, conda_env)
                except RuntimeError as e:
                    print(e)
                    print(f"Skipping combination: feature set '{feature_set}', tool '{tool}', database '{db}'")
