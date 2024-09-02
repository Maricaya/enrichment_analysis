import os
import subprocess
import yaml
import pandas as pd
import gseapy as gp

# Load the configuration file
config_path = os.path.abspath('test/config/example_enrichment_analysis_config.yaml')
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

result_path = os.path.abspath(config['result_path'])

# Load annotation file
annotation_path = os.path.abspath(config['annotation'])
annot = pd.read_csv(annotation_path, index_col='name')

# Extract gene sets and databases information
genes = annot.loc[annot['features_path'].str.endswith('.txt'), :]
genes_dict = genes.to_dict('index')

database_dict = config["local_databases"]
database_dict = {k: v for k, v in database_dict.items() if v != ""}

# Function to create the Conda environment if it doesn't exist
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

# Define the function to get gene set paths
def get_gene_path(gene_set):
    """
    Get the path to the gene set file for the given gene set.

    Parameters:
    - gene_set: The name of the gene set (str).

    Returns:
    - str: The absolute path to the gene set file.

    Raises:
    - ValueError: If the gene set is not found in the available dictionaries.
    """
    if gene_set in genes_dict.keys():
        return os.path.abspath(genes_dict[gene_set]['features_path'])
    else:
        raise ValueError(f"Gene set '{gene_set}' not found.")

# Define the main function to run ORA analysis using GSEApy
def run_ora_analysis(gene_set, database, conda_env):
    # Get the paths for the gene set and database
    gene_path = get_gene_path(gene_set)
    database_path = os.path.abspath(database_dict[database])

    output_dir = os.path.abspath(os.path.join(config['result_path'], 'enrichment_analysis', gene_set, 'ORA_GSEApy', database))
    result_file = os.path.join(output_dir, f'{gene_set}_{database}.csv')
    log_file = os.path.join('logs', f'gene_ORA_GSEApy_{gene_set}_{database}.log')

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Run ORA analysis using GSEApy
    try:
        # Perform Over-Representation Analysis (ORA) using GSEApy
        enrichment_results = gp.enrichr(
            gene_list=gene_path,
            gene_sets=database_path,
            outdir=output_dir
        )

        # Save the results to a CSV file
        enrichment_results.res2d.to_csv(result_file)

        print(f"ORA analysis completed successfully for gene set '{gene_set}' and database '{database}'.")

    except Exception as e:
        raise RuntimeError(f"ORA analysis failed for gene set '{gene_set}' and database '{database}'. Error: {e}")

if __name__ == "__main__":
    # Create Conda environment if it doesn't exist
    conda_env = 'gene_enrichment_analysis'
    env_file = os.path.abspath('workflow/envs/gene_enrichment_analysis.yaml')
    create_conda_env(conda_env, env_file)

    print(genes_dict.keys())
    print(database_dict.keys())

    # Example gene sets and databases
    for gene_set in genes_dict.keys():
        for database in database_dict.keys():
            print(f"Processing gene set: {gene_set} and database: {database}")
            run_ora_analysis(gene_set, database, conda_env)
