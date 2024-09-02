import os
import subprocess
import yaml

# Get the absolute path of the current script
current_dir = os.path.dirname(os.path.abspath(__file__))

# Load the configuration file
config_path = os.path.abspath('test/config/example_enrichment_analysis_config.yaml')
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

# Define the database dictionary
database_dict = {k: v for k, v in config["local_databases"].items() if v != ""}

# Function to get the user-provided local database path
def get_db_path(database):
    try:
        return database_dict[database]
    except KeyError:
        raise ValueError(f"Database '{database}' not found in the configuration.")

# Function to create the Conda environment if it doesn't exist
def create_conda_env(env_name, env_file):
    try:
        # Check if the environment already exists
        result = subprocess.run(['conda', 'env', 'list'], capture_output=True, text=True)
        if env_name in result.stdout:
            print(f"Conda environment '{env_name}' already exists.")
        else:
            # Create the environment
            print(f"Creating Conda environment '{env_name}'...")
            subprocess.run(['conda', 'env', 'create', '-f', env_file], check=True)
            print(f"Conda environment '{env_name}' created successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while creating the Conda environment: {e}")
        exit(1)

# Main function to prepare databases
def prepare_database(database):
    db_path = get_db_path(database)
    output_file = os.path.join("resources", config["project_name"], f"{database}.gmt")
    log_file = os.path.join("logs", "rules", f"prepare_databases_{database}.log")
    partition = config.get("partition")
    threads = config.get("threads", 1)
    mem_mb = config.get("mem", "16000")

    # Create the output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Define the command to run the script
    script_path = os.path.abspath("workflow/scripts/prepare_databases_GSEApy.py")
    if not os.path.exists(script_path):
        raise FileNotFoundError(f"Script {script_path} not found.")

    command = [
        "conda", "run", "-n", "gene_enrichment_analysis",  # Ensure the correct environment is activated
        "python", script_path,
        "--input", db_path,
        "--output", output_file,
        "--database", database  # Ensure the --database argument is included
    ]

    # Run the command and log the output
    with open(log_file, 'w') as log:
        result = subprocess.run(command, stdout=log, stderr=log)

    # Check if the command was successful
    if result.returncode != 0:
        print(f"An error occurred. Check the log file: {log_file}")
    else:
        print(f"Database preparation completed successfully. Output saved to: {output_file}")

# Prepare all databases defined in the configuration
if __name__ == "__main__":
    # Create the Conda environment if it doesn't exist
    env_file = os.path.abspath('workflow/envs/gene_enrichment_analysis.yaml')
    create_conda_env('gene_enrichment_analysis', env_file)

    for database in database_dict.keys():
        print(f"Preparing database: {database}")
        prepare_database(database)