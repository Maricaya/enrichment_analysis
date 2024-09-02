import os
import subprocess
import yaml
import pandas as pd

# 配置路径
result_path = "test/results"  # 示例中的结果路径，你可以根据实际路径进行修改

config_path = os.path.abspath('test/config/example_enrichment_analysis_config.yaml')
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)
# 生成绝对路径
result_path = os.path.abspath(result_path)

annotation_path = os.path.abspath(config['annotation'])
annot = pd.read_csv(annotation_path, index_col='name')
regions = annot.loc[annot['features_path'].str.endswith('.bed'), :]


regions_dict = regions.to_dict('index')

# 创建 Conda 环境（如果不存在）
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

# 处理 pycisTarget 结果的函数
def process_results_pycisTarget(region_set, database, conda_env):
    motif_hdf5_path = os.path.abspath(os.path.join(result_path, config['project_name'], region_set, 'pycisTarget', database, f'motif_enrichment_cistarget_{region_set}.hdf5'))
    motif_csv_path = os.path.abspath(os.path.join(result_path, config['project_name'], region_set, 'pycisTarget', database, f'{region_set}_{database}.csv'))

    print(f"Processing results for region set: {region_set} and database: {database}")

    # 检查输入文件是否存在
    if not os.path.exists(motif_hdf5_path):
        raise FileNotFoundError(f"Input file '{motif_hdf5_path}' not found.")

    try:
        # 使用 conda 运行 Python 脚本处理结果
        subprocess.run([
            'conda', 'run', '--no-capture-output', '--name', conda_env, 'python',
            './workflow/scripts/process_results_pycisTarget.py',
            motif_hdf5_path, motif_csv_path,
            'description_column',  # 这个参数需要根据你的实际需求修改
            region_set
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while processing the results: {e}")
        exit(1)

if __name__ == "__main__":
    # 假设你有一组区域集和数据库需要处理
    # region_sets = ["setCComplete", "setA", "setB"]  # 根据实际情况调整
    region_sets = regions_dict.keys()  # 从 regions_dict 中获取所有区域集
    print("region_sets", region_sets)
    # databases = ["hg38_screen_v10clust"]  # 根据实际情况调整
    # databases = config["local_databases"].keys()
    # 获取 pycistarget 数据库路径
    pycistarget_db_dict = config['pycistarget_parameters']['databases']
    databases = pycistarget_db_dict.keys()
    print("databases", databases)
    conda_env = "pycisTarget"

    # 创建 Conda 环境（如果不存在）
    env_file = os.path.abspath('workflow/envs/pycisTarget.yaml')
    create_conda_env(conda_env, env_file)

    # 为每个区域集和数据库运行分析
    for region_set in region_sets:
        for database in databases:
            process_results_pycisTarget(region_set, database, conda_env)
