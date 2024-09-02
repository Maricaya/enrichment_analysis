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

# 获取 pycistarget 数据库路径
pycistarget_db_dict = config['pycistarget_parameters']['databases']

# 定义辅助函数获取路径
def get_region_path(region_set):
    if region_set in regions_dict.keys():
        return os.path.abspath(regions_dict[region_set]['features_path'])
    else:
        raise ValueError(f"Region set '{region_set}' not found.")

def get_pycistarget_db_path(database):
    if database in pycistarget_db_dict.keys():
        return os.path.abspath(pycistarget_db_dict[database])
    else:
        raise ValueError(f"pycisTarget database '{database}' not found.")

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

# 执行 TFBS 基序富集分析的主函数
def region_motif_enrichment_analysis_pycisTarget(region_set, database):
    regions_path = get_region_path(region_set)
    ctx_db_path = get_pycistarget_db_path(database)
    motif2tf_path = os.path.abspath(config['pycistarget_parameters']['path_to_motif_annotations'])
    output_dir = os.path.abspath(os.path.join(config['result_path'], config['project_name'], region_set, 'pycisTarget', database))
    motif_hdf5 = os.path.join(output_dir, f"motif_enrichment_cistarget_{region_set}.hdf5")
    motif_html = os.path.join(output_dir, f"motif_enrichment_cistarget_{region_set}.html")
    log_file = os.path.abspath(os.path.join('logs', 'rules', f"region_enrichment_analysis_pycisTarget_{region_set}_{database}.log"))

    # 如果输出目录不存在，则创建
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    fraction_overlap = config['pycistarget_parameters']['fraction_overlap_w_cistarget_database']
    auc_threshold = config['pycistarget_parameters']['auc_threshold']
    nes_threshold = config['pycistarget_parameters']['nes_threshold']
    rank_threshold = config['pycistarget_parameters']['rank_threshold']
    annotation_version = config['pycistarget_parameters']['annotation_version']
    annotations_to_use = ",".join(config['pycistarget_parameters']['annotations_to_use'])
    motif_similarity_fdr = config['pycistarget_parameters']['motif_similarity_fdr']
    orthologous_identity_threshold = config['pycistarget_parameters']['orthologous_identity_threshold']
    species = 'homo_sapiens' if config['genome'] in ['hg19', 'hg38'] else 'mus_musculus' if config['genome'] in ['mm9', 'mm11'] else None

    # 构建命令
    command = (
        f"conda run -n pycisTarget "
        f"pycistarget cistarget "
        f"--cistarget_db_fname {ctx_db_path} "
        f"--bed_fname {regions_path} "
        f"--output_folder {output_dir} "
        f"--fr_overlap_w_ctx_db {fraction_overlap} "
        f"--auc_threshold {auc_threshold} "
        f"--nes_threshold {nes_threshold} "
        f"--rank_threshold {rank_threshold} "
        f"--path_to_motif_annotations {motif2tf_path} "
        f"--annotation_version {annotation_version} "
        f"--annotations_to_use {annotations_to_use} "
        f"--motif_similarity_fdr {motif_similarity_fdr} "
        f"--orthologous_identity_threshold {orthologous_identity_threshold} "
        f"--species {species} "
        f"--name {region_set} "
        f"--output_mode 'hdf5' "
        f"--write_html"
    )

    # 运行命令并记录输出
    with open(log_file, 'w') as log:
        result = subprocess.run(command, stdout=log, stderr=log, shell=True)

    # 检查命令是否成功执行
    if result.returncode != 0:
        print(f"An error occurred during the region TFBS motif enrichment analysis using pycisTarget. Check the log file: {log_file}")
        with open(motif_hdf5, 'w'), open(motif_html, 'w'):
            pass  # 创建空的输出文件以表示失败
    else:
        print(f"Region TFBS motif enrichment analysis completed successfully. Results saved to: {output_dir}")

# 为每个区域集和数据库运行分析
if __name__ == "__main__":
    # 创建 Conda 环境（如果不存在）
    env_file = os.path.abspath('workflow/envs/pycisTarget.yaml')
    create_conda_env('pycisTarget', env_file)

    # 假设有一组区域集合和数据库需要分析
    region_sets = regions_dict.keys()  # 从 regions_dict 中获取所有区域集
    databases = pycistarget_db_dict.keys()  # 从 pycistarget_db_dict 中获取所有数据库
    print("region_sets", region_sets)
    print("databases", databases)

    for region_set in region_sets:
        for database in databases:
            print(f"Running region motif enrichment analysis for region set: {region_set} and database: {database}")
            region_motif_enrichment_analysis_pycisTarget(region_set, database)
