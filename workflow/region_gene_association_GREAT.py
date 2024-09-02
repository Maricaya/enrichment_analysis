import os
import subprocess
import yaml
import pandas as pd

# 加载配置文件
config_path = os.path.abspath('test/config/example_enrichment_analysis_config.yaml')
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

# 加载 annotation 文件
annotation_path = os.path.abspath(config['annotation'])
annot = pd.read_csv(annotation_path, index_col='name')

# 提取 regions 信息
regions = annot.loc[annot['features_path'].str.endswith('.bed'), :]
regions_dict = regions.to_dict('index')

# 定义辅助函数获取路径
def get_region_path(region_set):
    if region_set in regions_dict.keys():
        return os.path.abspath(regions_dict[region_set]['features_path'])
    else:
        raise ValueError(f"Region set '{region_set}' not found.")

def get_first_database():
    if len(config["local_databases"]) > 0:
        return list(config["local_databases"].keys())[0]
    else:
        raise ValueError("No databases found in local_databases.")

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

# 执行 GREAT 基因关联分析的主函数
def region_gene_association_GREAT(region_set):
    regions_path = get_region_path(region_set)
    database = get_first_database()
    database_path = os.path.abspath(os.path.join("resources", config["project_name"], f"{database}.gmt"))
    output_dir = os.path.abspath(os.path.join(config['result_path'], config['project_name'], region_set, 'GREAT'))
    genes_output = os.path.join(output_dir, 'genes.txt')
    associations_table = os.path.join(output_dir, 'region_gene_associations.csv')
    associations_plot = os.path.join(output_dir, 'region_gene_associations.pdf')
    log_file = os.path.abspath(os.path.join('logs', 'rules', f"region_gene_association_GREAT_{region_set}.log"))

    # 如果输出目录不存在，则创建
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # 定义要运行的脚本命令
    script_path = os.path.abspath("workflow/scripts/region_gene_association_GREAT.R")
    if not os.path.exists(script_path):
        raise FileNotFoundError(f"Script {script_path} not found.")

    # 构建命令，传递参数
    command = (
        f"conda run -n region_enrichment_analysis Rscript {script_path} "
        f"{regions_path} "
        f"{database_path} "
        f"{genes_output} "
        f"{associations_table} "
        f"{associations_plot}"
    )

    # 运行命令并记录输出
    with open(log_file, 'w') as log:
        result = subprocess.run(command, stdout=log, stderr=log, shell=True)

    # 检查命令是否成功执行
    if result.returncode != 0:
        print(f"An error occurred during the region-gene association using GREAT. Check the log file: {log_file}")
    else:
        print(f"Region-gene association completed successfully. Results saved to: {output_dir}")

# 为每个区域集运行分析
if __name__ == "__main__":
    # 创建 Conda 环境（如果不存在）
    env_file = os.path.abspath('workflow/envs/region_enrichment_analysis.yaml')
    create_conda_env('region_enrichment_analysis', env_file)

    # 假设有一组区域集合需要分析
    region_sets = regions_dict.keys()  # 从 regions_dict 中获取所有区域集

    for region_set in region_sets:
        print(f"Running region-gene association for region set: {region_set}")
        region_gene_association_GREAT(region_set)
