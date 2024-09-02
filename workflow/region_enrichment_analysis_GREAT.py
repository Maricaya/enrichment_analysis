import os
import subprocess
import yaml
import pandas as pd

# Load the configuration file
config_path = os.path.abspath('test/config/example_enrichment_analysis_config.yaml')
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

result_path = os.path.abspath(config['result_path'])

# Load annotation file
annotation_path = os.path.abspath(config['annotation'])
annot = pd.read_csv(annotation_path, index_col='name')

# Extract regions and databases information
regions = annot.loc[annot['features_path'].str.endswith('.bed'), :]
regions_dict = regions.to_dict('index')

background_regions = annot.loc[:, ['background_name', 'background_path']]
background_regions_dict = background_regions.drop_duplicates().set_index('background_name').to_dict('index')

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

# Define the function to get region paths
def get_region_path(region_set):
    """
    Get the path to the region file for the given region set.

    Parameters:
    - region_set: The name of the region set (str).

    Returns:
    - str: The absolute path to the region file.

    Raises:
    - ValueError: If the region set is not found in the available dictionaries.
    """
    if region_set in regions_dict.keys():
        return os.path.abspath(regions_dict[region_set]['features_path'])
    elif region_set in background_regions_dict.keys():
        return os.path.abspath(background_regions_dict[region_set]['background_path'])
    else:
        raise ValueError(f"Region set '{region_set}' not found.")

def run_great_analysis(region_set, database, conda_env):
    # Get the paths for the region set and database
    region_path = get_region_path(region_set)
    background_path = get_region_path(region_set)
    database_path = os.path.abspath(os.path.join("resources", config["project_name"], f"{database}.gmt"))

    output_dir = os.path.abspath(os.path.join(config['result_path'], 'enrichment_analysis', region_set, 'GREAT', database))
    result_file = os.path.join(output_dir, f'{region_set}_{database}.csv')
    log_file = os.path.join('logs', f'region_enrichment_analysis_GREAT_{region_set}_{database}.log')

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Command to run GREAT analysis within the specified Conda environment
    command = [
        'conda', 'run', '--no-capture-output', '--name', conda_env,
        'Rscript', 'workflow/scripts/region_enrichment_analysis_GREAT.R',
        region_path, background_path, database_path, result_file, config['genome'], str(config.get('threads', 1))
    ]

    # Run the command
    with open(log_file, 'w') as log:
        result = subprocess.run(command, stdout=log, stderr=log, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"GREAT analysis failed for region set '{region_set}' and database '{database}'. Check the log file {log_file} for details.")
    else:
        print(f"GREAT analysis completed successfully for region set '{region_set}' and database '{database}'.")

if __name__ == "__main__":
    # Create Conda environment if it doesn't exist
    conda_env = 'region_enrichment_analysis'
    env_file = os.path.abspath('workflow/envs/region_enrichment_analysis.yaml')
    create_conda_env(conda_env, env_file)

    # Example region sets and databases
    for region_set in regions_dict.keys():
        for database in config['local_databases'].keys():
            print(f"Processing region set: {region_set} and database: {database}")
            run_great_analysis(region_set, database, conda_env)
