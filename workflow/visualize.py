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

# 自动获取 group、tool 和 db 的函数
def get_groups_tools_dbs(config):
    annotation_path = os.path.abspath(config['annotation'])
    annot = pd.read_csv(annotation_path, index_col='name')

    groups = annot['group'].unique().tolist()
    tools = list(config["column_names"].keys())
    databases = list(config["local_databases"].keys())

    return groups, tools, databases

# 主函数
def main():
    conda_env = 'visualization'
    if conda_env is None:
        raise RuntimeError("Conda environment not found. Ensure the script is run within a Snakemake conda environment.")

    groups, tools, databases = get_groups_tools_dbs(config)

    for group in groups:
        for tool in tools:
            for db in databases:
                results_all = os.path.join(result_path, "enrichment_analysis", group, tool, db, f"{group}_{db}_all.csv")
                summary_plot = os.path.join(result_path, "enrichment_analysis", group, tool, db, f"{group}_{db}_summary.png")
                adjp_hm = os.path.join(result_path, "enrichment_analysis", group, tool, db, f"{group}_{db}_adjp_heatmap.pdf")
                effect_hm = os.path.join(result_path, "enrichment_analysis", group, tool, db, f"{group}_{db}_effect_heatmap.pdf")
                log_file = os.path.join("logs", f"visualize_{group}_{tool}_{db}.log")

                # 创建日志和输出目录
                os.makedirs(os.path.dirname(summary_plot), exist_ok=True)
                os.makedirs(os.path.dirname(adjp_hm), exist_ok=True)
                os.makedirs(os.path.dirname(effect_hm), exist_ok=True)
                os.makedirs(os.path.dirname(log_file), exist_ok=True)

                # 调用 R 脚本生成可视化
                script_path = os.path.abspath("workflow/scripts/overview_plot.R")
                command = [
                    'conda', 'run', '--no-capture-output', '--name', conda_env,
                    'Rscript', script_path,
                    results_all, summary_plot, adjp_hm, effect_hm,
                    tool, db, group, config_path
                ]

                with open(log_file, 'w') as log:
                    result = subprocess.run(command, stdout=log, stderr=log, text=True)

                # 检查执行结果
                if result.returncode != 0:
                    print(f"Visualization failed for group '{group}', tool '{tool}', and database '{db}'. Check the log file {log_file} for details.")
                else:
                    print(f"Visualization completed successfully for group '{group}', tool '{tool}', and database '{db}'.")
                    print(f"Summary plot saved in {summary_plot}")
                    print(f"AdjP heatmap saved in {adjp_hm}")
                    print(f"Effect heatmap saved in {effect_hm}")

if __name__ == "__main__":
    main()
