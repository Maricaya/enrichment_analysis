def process_results_pycisTarget(region_set, database, term_col, conda_env):
    result_path = os.path.abspath('test/results')

    # Generate paths
    motif_hdf5_path = os.path.join(result_path, region_set, 'pycisTarget', database, f'motif_enrichment_cistarget_{region_set}.hdf5')
    motif_csv_path = os.path.join(result_path, region_set, 'pycisTarget', database, f'{region_set}_{database}.csv')

    print(f"Processing results for region set: {region_set} and database: {database}")
    print(f"Expected .hdf5 file path: {motif_hdf5_path}")

    if not os.path.exists(motif_hdf5_path):
        print(f"File not found at: {motif_hdf5_path}")
        raise FileNotFoundError(f"Input file '{motif_hdf5_path}' not found.")

    # Load results from HDF5
    results = read_hdf5(motif_hdf5_path)
    results_df = results[region_set].motif_enrichment

    # Print available columns for debugging
    print("Available columns:", results_df.columns)

    # Ensure term_col exists in DataFrame
    if term_col not in results_df.columns:
        raise KeyError(f"Column '{term_col}' not found in the DataFrame.")

    # Add 'description' column to DataFrame
    results_df["description"] = results_df["motif"] + "(" + results_df[term_col] + ")"

    # Save to CSV
    results_df.to_csv(motif_csv_path, index=False)
