import os
import subprocess
import yaml
import pandas as pd
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr

# Enable the conversion between R and pandas data frames
pandas2ri.activate()

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

rcistarget_db_dict = config["rcistarget_parameters"]["databases"]
rcistarget_db_dict = {k: v for k, v in rcistarget_db_dict.items() if v != ""}

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

# Define the main function to run RcisTarget analysis
def run_rcistarget_analysis(gene_set, database, conda_env):
    # Get the paths for the gene set and database
    gene_path = get_gene_path(gene_set)
    database_path = os.path.abspath(rcistarget_db_dict[database])
    motif_annotation = os.path.abspath(config["rcistarget_parameters"]["motifAnnot"])

    output_dir = os.path.abspath(os.path.join(config['result_path'], 'enrichment_analysis', gene_set, 'RcisTarget', database))
    result_file = os.path.join(output_dir, f'{gene_set}_{database}.csv')
    log_file = os.path.join('logs', f'gene_motif_enrichment_analysis_RcisTarget_{gene_set}_{database}.log')

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Command to run RcisTarget analysis within the specified Conda environment
    command = [
        'conda', 'run', '--no-capture-output', '--name', conda_env,
        'Rscript', 'workflow/scripts/gene_enrichment_analysis_RcisTarget.R',
        gene_path, database_path, motif_annotation, result_file
    ]

    # Run the command
    with open(log_file, 'w') as log:
        result = subprocess.run(command, stdout=log, stderr=log, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"RcisTarget analysis failed for gene set '{gene_set}' and database '{database}'. Check the log file {log_file} for details.")
    else:
        print(f"RcisTarget analysis completed successfully for gene set '{gene_set}' and database '{database}'.")

if __name__ == "__main__":
    # Create Conda environment if it doesn't exist
    conda_env = 'RcisTarget'
    env_file = os.path.abspath('workflow/envs/RcisTarget.yaml')
    create_conda_env(conda_env, env_file)

    print(genes_dict.keys())
    print(rcistarget_db_dict.keys())

    # Example gene sets and databases
    for gene_set in genes_dict.keys():
        for database in rcistarget_db_dict.keys():
            print(f"Processing gene set: {gene_set} and database: {database}")
            run_rcistarget_analysis(gene_set, database, conda_env)
