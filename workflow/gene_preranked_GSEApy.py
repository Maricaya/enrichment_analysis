import os
import pandas as pd
import gseapy as gp
import yaml

# 加载配置文件
config_path = os.path.abspath('test/config/example_enrichment_analysis_config.yaml')
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

# 加载 annotation 文件
annotation_path = os.path.abspath(config['annotation'])

print("annotation_path", annotation_path)
annot = pd.read_csv(annotation_path, index_col='name')

print("annot", annot)

# 提取 genes 和 regions 信息
genes = annot.loc[annot['features_path'].str.endswith('.txt'), :]
print("genes", genes)
regions = annot.loc[annot['features_path'].str.endswith('.bed'), :]

genes_dict = genes.to_dict('index')
regions_dict = regions.to_dict('index')
result_path = os.path.abspath(config['result_path'])

# 定义获取基因集路径的函数
def get_gene_path(gene_set):
    if gene_set in genes_dict:
        return os.path.abspath(genes_dict[gene_set]['features_path'])
    elif gene_set in regions_dict:
        return os.path.abspath(os.path.join(result_path, gene_set, 'GREAT', 'genes.txt'))
    else:
        raise ValueError(f"Gene set '{gene_set}' not found.")

# 定义获取数据库路径的函数
def get_database():
    if len(config["local_databases"]) > 0:
        return list(config["local_databases"].keys())
    else:
        raise ValueError("No databases found in local_databases.")

# 定义进行 GSEA 分析的主函数
def gene_preranked_GSEApy(gene_set, database):
    ranked_genes_path = get_gene_path(gene_set)
    database_path = os.path.abspath(os.path.join("resources", config["project_name"], f"{database}.gmt"))
    output_dir = os.path.abspath(os.path.join(result_path, gene_set, 'preranked_GSEApy', database))
    report_dir = os.path.join(output_dir, 'report')
    log_file = os.path.abspath(os.path.join('logs', 'rules', f"gene_preranked_GSEApy_{gene_set}_{database}.log"))

    # 如果输出目录不存在，则创建
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # 加载基因排名数据
    preranked_data = pd.read_csv(ranked_genes_path, sep="\t", index_col=0)

    # 运行 GSEApy preranked
    gsea_results = gp.prerank(
        rnk=preranked_data,
        gene_sets=database_path,
        outdir=report_dir,
        format='png',
        seed=6,
        verbose=True
    )

    # 保存重要的 GSEApy 结果
    gsea_results.res2d.to_csv(os.path.join(output_dir, f"{gene_set}_{database}_GSEApy_results.csv"))

    # 绘制前几个结果的富集图
    for i, term in enumerate(gsea_results.res2d.index[:5]):
        gseaplot(gsea_results.ranking, term=term, ofname=os.path.join(output_dir, f"gseaplot_{i}.png"))

# 为每个基因集运行分析
if __name__ == "__main__":
    gene_sets = genes_dict.keys()  # 从 genes_dict 中获取所有基因集
    databases = get_database()  # 获取所有数据库

    print("gene_sets", gene_sets)
    print("databases", databases)

    for gene_set in gene_sets:
        for database in databases:
            print(f"Running GSEA for gene set: {gene_set} and database: {database}")
            gene_preranked_GSEApy(gene_set, database)
