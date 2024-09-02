#!/usr/bin/env python3

import os
import subprocess
import yaml
import pandas as pd

# 加载配置文件
config_path = os.path.abspath('test/config/example_enrichment_analysis_config.yaml')
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

result_path = os.path.abspath(config['result_path'])

# 定义从组获取路径的函数
def get_group_paths(group, tool, db, annot, regions_dict, genes_dict, rnk_dict):
    feature_sets = list(annot.index[annot["group"] == group])

    # 根据工具选择特定的 feature_sets
    if tool in ["GREAT", "LOLA", "pycisTarget"]:
        feature_sets = [fs for fs in feature_sets if fs in regions_dict.keys()]
    elif tool in ["ORA_GSEApy", "RcisTarget"]:
        feature_sets = [fs for fs in feature_sets if fs in genes_dict.keys()] + \
                       [fs for fs in feature_sets if fs in regions_dict.keys()]
    elif tool == "preranked_GSEApy":
        feature_sets = [fs for fs in feature_sets if fs in rnk_dict.keys()]

    # 生成文件路径
    paths = [
        os.path.join(result_path, fs, tool, db, f"{fs}_{db}.csv")
        for fs in feature_sets
    ]
    return paths

# 自动获取 group、tool 和 db 的函数
def get_groups_tools_dbs(config):
    annotation_path = os.path.abspath(config['annotation'])
    annot = pd.read_csv(annotation_path, index_col='name')

    groups = annot['group'].unique().tolist()
    tools = list(config["column_names"].keys())
    databases = list(config["local_databases"].keys())

    return groups, tools, databases, annot

# 从配置文件加载字典
def load_dictionaries(config):
    annot = pd.read_csv(config["annotation"], index_col='name')

    genes = annot.loc[annot['features_path'].str.endswith('.txt'), :]
    regions = annot.loc[annot['features_path'].str.endswith('.bed'), :]
    genes_dict = genes.to_dict('index')
    regions_dict = regions.to_dict('index')

    rnk = annot.loc[annot['features_path'].str.endswith('.csv'),:]
    rnk_dict = rnk.to_dict('index')

    return regions_dict, genes_dict, rnk_dict

# 主函数
def main():
    conda_env = os.environ.get('CONDA_PREFIX')
    if conda_env is None:
        raise RuntimeError("Conda environment not found. Ensure the script is run within a Snakemake conda environment.")

    # 加载注释文件和字典
    groups, tools, databases, annot = get_groups_tools_dbs(config)
    regions_dict, genes_dict, rnk_dict = load_dictionaries(config)

    for group in groups:
        for tool in tools:
            for db in databases:
                enrichment_results = get_group_paths(group, tool, db, annot, regions_dict, genes_dict, rnk_dict)
                if not enrichment_results:
                    continue

                results_all = os.path.join(result_path, "enrichment_analysis", group, tool, db, f"{group}_{db}_all.csv")
                results_sig = os.path.join(result_path, "enrichment_analysis", group, tool, db, f"{group}_{db}_sig.csv")
                log_file = os.path.join("logs", f"aggregate_{group}_{tool}_{db}.log")

                os.makedirs(os.path.dirname(results_all), exist_ok=True)
                os.makedirs(os.path.dirname(log_file), exist_ok=True)

                # 调用 scripts/aggregate.py
                script_path = os.path.abspath("workflow/scripts/aggregate.py")
                command = [
                    'conda', 'run', '--no-capture-output', '-p', conda_env,
                    'python', script_path,
                    '--enrichment_results', *enrichment_results,
                    '--results_all', results_all,
                    '--results_sig', results_sig,
                    '--group', group,
                    '--tool', tool,
                    '--db', db,
                    '--config', config_path
                ]

                with open(log_file, 'w') as log:
                    result = subprocess.run(command, stdout=log, stderr=log, text=True)

                if result.returncode != 0:
                    print(f"Aggregation failed for group '{group}', tool '{tool}', and database '{db}'. Check the log file {log_file} for details.")
                else:
                    print(f"Aggregation completed successfully for group '{group}', tool '{tool}', and database '{db}'.")
                    print(f"Results saved in:\n  - {results_all}\n  - {results_sig}")

if __name__ == "__main__":
    main()
