#!/usr/bin/env python3

import os
import sys
import yaml
import pandas as pd
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description="Aggregate enrichment results.")
parser.add_argument('--enrichment_results', nargs='+', required=True, help="List of enrichment result files.")
parser.add_argument('--results_all', required=True, help="Path to save all combined results.")
parser.add_argument('--results_sig', required=True, help="Path to save significant results.")
parser.add_argument('--group', required=True, help="Group name.")
parser.add_argument('--tool', required=True, help="Tool name.")
parser.add_argument('--db', required=True, help="Database name.")
parser.add_argument('--config', required=True, help="Path to the config file.")
args = parser.parse_args()

# 加载配置文件
with open(args.config, 'r') as file:
    config_data = yaml.safe_load(file)

term_col = config_data["column_names"][args.tool]["term"]
adjp_col = config_data["column_names"][args.tool]["adj_pvalue"]
adjp_th = config_data["adjp_th"][args.tool]

# 加载所有的结果文件
results_list = []
for result_path in args.enrichment_results:
    if os.path.exists(result_path) and os.path.getsize(result_path) > 0:
        tmp_name = os.path.basename(result_path).replace(f"_{args.db}.csv", "")
        tmp_res = pd.read_csv(result_path, index_col=0)
        tmp_res['name'] = tmp_name
        results_list.append(tmp_res)

# 如果没有有效的结果文件，创建空文件并退出
if not results_list:
    pd.DataFrame().to_csv(args.results_all)
    pd.DataFrame().to_csv(args.results_sig)
    sys.exit(0)

# 将所有结果文件合并为一个 DataFrame
result_df = pd.concat(results_list, axis=0)
result_df.to_csv(args.results_all)  # 保存所有的合并结果

# 根据显著性水平过滤结果
if args.tool in ["pycisTarget", "RcisTarget"]:
    sig_terms = result_df.loc[result_df[adjp_col] >= adjp_th, term_col].unique()
else:
    sig_terms = result_df.loc[result_df[adjp_col] <= adjp_th, term_col].unique()

result_sig_df = result_df.loc[result_df[term_col].isin(sig_terms), :]
result_sig_df.to_csv(args.results_sig)  # 保存显著性过滤后的结果
