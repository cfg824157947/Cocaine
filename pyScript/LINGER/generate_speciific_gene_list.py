import pandas as pd
# Assuming `df` is your DataFrame
# Filter terms where value is smaller than 0.05 for each cell type


def filter_terms(df):
    filtered_terms = {}
    for cell_type in df.columns:
        filtered = df[df[cell_type]]
        filtered_terms[cell_type] = filtered.index.tolist()
    max_length = max(len(terms) for terms in filtered_terms.values())  # Find the maximum list length
    filtered_terms_df = pd.DataFrame({
        cell_type: terms + [None] * (max_length - len(terms))  # Pad shorter lists with None
        for cell_type, terms in filtered_terms.items()
    })
    return filtered_terms_df




file_list = [item for item in os.listdir(driver_result_dir) if ('C_result' in item)]
for entry in file_list:
    driver_result_file = os.path.join(driver_result_dir, entry)
    print(driver_result_file)
    Cdf = pd.read_csv(driver_result_file, index_col=0)
    Qdf = pd.read_csv(driver_result_file.replace("C_result", "Q_result"), index_col=0)
    up_down_list = ['up', 'down']
    for up_down in up_down_list:
        if up_down == 'up':
            selected_elements_df = ((Cdf > 0.1) & (Qdf < 0.01))
        else:
            selected_elements_df = ((Cdf < -0.1) & (Qdf < 0.01))
        filtered_terms_df = filter_terms(selected_elements_df)
        # Save the DataFrame to a CSV file
        result_file = driver_result_file.replace("C_result", ("filtered_terms_"+up_down))
        filtered_terms_df.to_csv(result_file, index=False)

# Save the DataFrame to a CSV file
filtered_terms_df.to_csv("filtered_GSEA_terms.csv", index=False)

print("Filtered terms have been saved to 'filtered_terms.csv'.")
