import sys
import os
import argparse
import pandas as pd
from scipy.stats import ttest_ind

# Setup argument parser
parser = argparse.ArgumentParser(description='miRNA Analyzer Script')
parser.add_argument('file_path', type=str, help='Path to the CSV file containing miRNA data')
parser.add_argument('--stablehk', type=str, help='Space-separated list of stable housekeeping miRNAs')
parser.add_argument('--stableneg', type=str, help='Space-separated list of stable negative miRNAs')
parser.add_argument('--stablepos', type=str, help='Space-separated list of stable positive miRNAs')

# Parse arguments
args = parser.parse_args()

# Function to prompt for input or use default
def prompt_for_mirnas(prompt, default):
    return input(f"{prompt} [default: {', '.join(default)}]: ") or default

# Default lists of miRNAs
default_stable_hk = ['RPL19', 'RPLP0', 'GAPDH', 'ACTB']
default_stable_neg = ['NEG_G', 'NEG_D', 'NEG_A', 'NEG_B']
default_stable_pos = ['POS_C', 'POS_F', 'POS_B', 'POS_A']

# File path is a required argument
file_path = args.file_path

# Set the identified "stable" reference miRNAs (found via external tools such as RefFinder ( http://blooge.cn/RefFinder/ ))
# Use arguments if provided, otherwise prompt
stable_housekeeping_miRNAs = args.stablehk.split() if args.stablehk else prompt_for_mirnas("Stable Housekeeping miRNAs", default_stable_hk)
stable_negative_miRNAs = args.stableneg.split() if args.stableneg else prompt_for_mirnas("Stable Negative miRNAs", default_stable_neg)
stable_positive_miRNAs = args.stablepos.split() if args.stablepos else prompt_for_mirnas("Stable Positive miRNAs", default_stable_pos)

# Load the data from the CSV file provided at the command line
data = pd.read_csv(file_path)

# Display the first few rows of the dataset to understand its structure
data.head()


# Identifying the stable housekeeping miRNAs in the dataset
all_housekeeping_miRNAs = data[data['Code Class'] == 'Housekeeping']
housekeeping_miRNAs = all_housekeeping_miRNAs[all_housekeeping_miRNAs['Name'].isin(stable_housekeeping_miRNAs)]

# Checking the first and last few rows of the housekeeping miRNAs to confirm
print(housekeeping_miRNAs)


# Identifying the stable negative miRNAs in the dataset
all_negative_miRNAs = data[data['Code Class'] == 'Negative']
negative_miRNAs = all_negative_miRNAs[all_negative_miRNAs['Name'].isin(stable_negative_miRNAs)]

# Checking the first and last few rows of the negative miRNAs to confirm
print(negative_miRNAs)


# Identifying the positive miRNAs in the dataset
all_positive_miRNAs = data[data['Code Class'] == 'Positive']
positive_miRNAs = all_positive_miRNAs[all_positive_miRNAs['Name'].isin(stable_positive_miRNAs)]

# Checking the first and last few rows of the positive miRNAs to confirm
print(positive_miRNAs)


# Exclude the housekeeping, negative, and positive to create a clean miRNA dataset
invalid_code_classes = ['Housekeeping', 'Negative', 'Positive']
clean_miRNAs = data[~data['Code Class'].isin(invalid_code_classes)] 

# Checking the first and last few rows of the clean miRNAs to confirm
print(clean_miRNAs)


# Extracting the expression data for housekeeping miRNAs
housekeeping_expression = housekeeping_miRNAs.iloc[:, 3:]

# Extracting the expression data for negative miRNAs
negative_expression = negative_miRNAs.iloc[:, 3:]

# Extracting the expression data for positive miRNAs
positive_expression = positive_miRNAs.iloc[:, 3:]

# Calculating the average expression of housekeeping miRNAs for each sample
average_housekeeping_expression = housekeeping_expression.mean()

# Calculating the average expression of negative miRNAs for each sample
average_negative_expression = negative_expression.mean()

# Calculating the average expression of negative miRNAs for each sample
average_positive_expression = positive_expression.mean()

# Normalizing the expression levels of all miRNAs by the average housekeeping expression
housekeeping_normalized_data = clean_miRNAs.iloc[:, 3:].div(average_housekeeping_expression)

# Normalizing the expression levels of all clean miRNAs by the average negative expression
negative_normalized_data = clean_miRNAs.iloc[:, 3:].div(average_negative_expression)

# Normalizing the expression levels of all clean miRNAs by the average positive expression
positive_normalized_data = clean_miRNAs.iloc[:, 3:].div(average_positive_expression)

# Adding back the non-expression columns to the housekeeping normalized dataset
housekeeping_normalized_data.insert(0, 'Accession', data['Accession'])
housekeeping_normalized_data.insert(0, 'Name', data['Name'])
housekeeping_normalized_data.insert(0, 'Code Class', data['Code Class'])

# Displaying the first few rows of the housekeeping normalized dataset
print(housekeeping_normalized_data)

# Adding back the non-expression columns to the negative normalized dataset
negative_normalized_data.insert(0, 'Accession', data['Accession'])
negative_normalized_data.insert(0, 'Name', data['Name'])
negative_normalized_data.insert(0, 'Code Class', data['Code Class'])

# Displaying the first few rows of the negative normalized dataset
print(negative_normalized_data)

# Adding back the non-expression columns to the positive normalized dataset
positive_normalized_data.insert(0, 'Accession', data['Accession'])
positive_normalized_data.insert(0, 'Name', data['Name'])
positive_normalized_data.insert(0, 'Code Class', data['Code Class'])

# Displaying the first few rows of the positive normalized dataset
print(positive_normalized_data)

# Extracting columns for each housekeeping group
hk_af_columns = [col for col in housekeeping_normalized_data.columns if col.startswith('AF')]
hk_am_columns = [col for col in housekeeping_normalized_data.columns if col.startswith('AM')]
hk_fc_columns = [col for col in housekeeping_normalized_data.columns if col.startswith('FC')]
hk_mc_columns = [col for col in housekeeping_normalized_data.columns if col.startswith('MC')]

# Extracting columns for each negative group
neg_af_columns = [col for col in negative_normalized_data.columns if col.startswith('AF')]
neg_am_columns = [col for col in negative_normalized_data.columns if col.startswith('AM')]
neg_fc_columns = [col for col in negative_normalized_data.columns if col.startswith('FC')]
neg_mc_columns = [col for col in negative_normalized_data.columns if col.startswith('MC')]

# Extracting columns for each positive group
pos_af_columns = [col for col in positive_normalized_data.columns if col.startswith('AF')]
pos_am_columns = [col for col in positive_normalized_data.columns if col.startswith('AM')]
pos_fc_columns = [col for col in positive_normalized_data.columns if col.startswith('FC')]
pos_mc_columns = [col for col in positive_normalized_data.columns if col.startswith('MC')]

# Calculating the mean for each group
grouped_data = pd.DataFrame({
    'Code Class': clean_miRNAs['Code Class'],
    'Name': clean_miRNAs['Name'],
    'Accession': clean_miRNAs['Accession'],
    'AF_Mean_HK': housekeeping_normalized_data[hk_af_columns].mean(axis=1),
    'AM_Mean_HK': housekeeping_normalized_data[hk_am_columns].mean(axis=1),
    'FC_Mean_HK': housekeeping_normalized_data[hk_fc_columns].mean(axis=1),
    'MC_Mean_HK': housekeeping_normalized_data[hk_mc_columns].mean(axis=1),
    'AF_Mean_NEG': negative_normalized_data[neg_af_columns].mean(axis=1),
    'AM_Mean_NEG': negative_normalized_data[neg_am_columns].mean(axis=1),
    'FC_Mean_NEG': negative_normalized_data[neg_fc_columns].mean(axis=1),
    'MC_Mean_NEG': negative_normalized_data[neg_mc_columns].mean(axis=1),
    'AF_Mean_POS': positive_normalized_data[pos_af_columns].mean(axis=1),
    'AM_Mean_POS': positive_normalized_data[pos_am_columns].mean(axis=1),
    'FC_Mean_POS': positive_normalized_data[pos_fc_columns].mean(axis=1),
    'MC_Mean_POS': positive_normalized_data[pos_mc_columns].mean(axis=1)
})

# Displaying the first few rows of the grouped data
print(grouped_data)

# Calculating fold changes for the specified comparisons
grouped_data['AF_to_AM_FoldChange_HK'] = grouped_data['AF_Mean_HK'] / grouped_data['AM_Mean_HK']
grouped_data['AF_to_FC_FoldChange_HK'] = grouped_data['AF_Mean_HK'] / grouped_data['FC_Mean_HK']
grouped_data['AM_to_MC_FoldChange_HK'] = grouped_data['AM_Mean_HK'] / grouped_data['MC_Mean_HK']
grouped_data['FC_to_MC_FoldChange_HK'] = grouped_data['FC_Mean_HK'] / grouped_data['MC_Mean_HK']
grouped_data['AF_to_AM_FoldChange_NEG'] = grouped_data['AF_Mean_NEG'] / grouped_data['AM_Mean_NEG']
grouped_data['AF_to_FC_FoldChange_NEG'] = grouped_data['AF_Mean_NEG'] / grouped_data['FC_Mean_NEG']
grouped_data['AM_to_MC_FoldChange_NEG'] = grouped_data['AM_Mean_NEG'] / grouped_data['MC_Mean_NEG']
grouped_data['FC_to_MC_FoldChange_NEG'] = grouped_data['FC_Mean_NEG'] / grouped_data['MC_Mean_NEG']
grouped_data['AF_to_AM_FoldChange_POS'] = grouped_data['AF_Mean_POS'] / grouped_data['AM_Mean_POS']
grouped_data['AF_to_FC_FoldChange_POS'] = grouped_data['AF_Mean_POS'] / grouped_data['FC_Mean_POS']
grouped_data['AM_to_MC_FoldChange_POS'] = grouped_data['AM_Mean_POS'] / grouped_data['MC_Mean_POS']
grouped_data['FC_to_MC_FoldChange_POS'] = grouped_data['FC_Mean_POS'] / grouped_data['MC_Mean_POS']

# Displaying the first few rows of the grouped data with the fold changes
print(grouped_data)


# Function to perform t-test between two normalized groups
def perform_ttest(data_normalized, group1, group2):
    p_values = []
    for index, row in data_normalized.iterrows():
        stat, p_val = ttest_ind(row[group1].astype(str).astype(float), row[group2].astype(str).astype(float), nan_policy='omit')
        p_values.append(p_val)
    return p_values

# Calculating p-values for each comparison
grouped_data['AF_AM_p_value_HK'] = perform_ttest(housekeeping_normalized_data, hk_af_columns, hk_am_columns)
grouped_data['AF_FC_p_value_HK'] = perform_ttest(housekeeping_normalized_data, hk_af_columns, hk_fc_columns)
grouped_data['AM_MC_p_value_HK'] = perform_ttest(housekeeping_normalized_data, hk_am_columns, hk_mc_columns)
grouped_data['FC_MC_p_value_HK'] = perform_ttest(housekeeping_normalized_data, hk_fc_columns, hk_mc_columns)
grouped_data['AF_AM_p_value_NEG'] = perform_ttest(negative_normalized_data, neg_af_columns, neg_am_columns)
grouped_data['AF_FC_p_value_NEG'] = perform_ttest(negative_normalized_data, neg_af_columns, neg_fc_columns)
grouped_data['AM_MC_p_value_NEG'] = perform_ttest(negative_normalized_data, neg_am_columns, neg_mc_columns)
grouped_data['FC_MC_p_value_NEG'] = perform_ttest(negative_normalized_data, neg_fc_columns, neg_mc_columns)
grouped_data['AF_AM_p_value_POS'] = perform_ttest(positive_normalized_data, pos_af_columns, pos_am_columns)
grouped_data['AF_FC_p_value_POS'] = perform_ttest(positive_normalized_data, pos_af_columns, pos_fc_columns)
grouped_data['AM_MC_p_value_POS'] = perform_ttest(positive_normalized_data, pos_am_columns, pos_mc_columns)
grouped_data['FC_MC_p_value_POS'] = perform_ttest(positive_normalized_data, pos_fc_columns, pos_mc_columns)

# Displaying the first few rows of the grouped data with the fold changes and p-values
print(grouped_data)

# Do some quick analysis to check the data for significant miRNAs
# Criteria for selection
fold_change_threshold = 2
p_value_threshold = 0.1

# Function to select miRNAs based on fold change and p-value criteria
def select_significant_miRNAs(data_grouped, fold_change_column, p_value_column):
    significant_miRNAs = data_grouped[
        (data_grouped[fold_change_column].abs() >= fold_change_threshold) & 
        (data_grouped[p_value_column] < p_value_threshold)
    ]
    return significant_miRNAs

# Selecting significant miRNAs for each comparison
significant_AF_AM_HK = select_significant_miRNAs(grouped_data, 'AF_to_AM_FoldChange_HK', 'AF_AM_p_value_HK')
significant_AF_FC_HK = select_significant_miRNAs(grouped_data, 'AF_to_FC_FoldChange_HK', 'AF_FC_p_value_HK')
significant_AM_MC_HK = select_significant_miRNAs(grouped_data, 'AM_to_MC_FoldChange_HK', 'AM_MC_p_value_HK')
significant_FC_MC_HK = select_significant_miRNAs(grouped_data, 'FC_to_MC_FoldChange_HK', 'FC_MC_p_value_HK')
significant_AF_AM_NEG = select_significant_miRNAs(grouped_data, 'AF_to_AM_FoldChange_NEG', 'AF_AM_p_value_NEG')
significant_AF_FC_NEG = select_significant_miRNAs(grouped_data, 'AF_to_FC_FoldChange_NEG', 'AF_FC_p_value_NEG')
significant_AM_MC_NEG = select_significant_miRNAs(grouped_data, 'AM_to_MC_FoldChange_NEG', 'AM_MC_p_value_NEG')
significant_FC_MC_NEG = select_significant_miRNAs(grouped_data, 'FC_to_MC_FoldChange_NEG', 'FC_MC_p_value_NEG')
significant_AF_AM_POS = select_significant_miRNAs(grouped_data, 'AF_to_AM_FoldChange_POS', 'AF_AM_p_value_POS')
significant_AF_FC_POS = select_significant_miRNAs(grouped_data, 'AF_to_FC_FoldChange_POS', 'AF_FC_p_value_POS')
significant_AM_MC_POS = select_significant_miRNAs(grouped_data, 'AM_to_MC_FoldChange_POS', 'AM_MC_p_value_POS')
significant_FC_MC_POS = select_significant_miRNAs(grouped_data, 'FC_to_MC_FoldChange_POS', 'FC_MC_p_value_POS')

# Summarizing the results
significant_summary = {
    'AF_vs_AM_HK': significant_AF_AM_HK[['Name', 'AF_to_AM_FoldChange_HK', 'AF_AM_p_value_HK']],
    'AF_vs_FC_HK': significant_AF_FC_HK[['Name', 'AF_to_FC_FoldChange_HK', 'AF_FC_p_value_HK']],
    'AM_vs_MC_HK': significant_AM_MC_HK[['Name', 'AM_to_MC_FoldChange_HK', 'AM_MC_p_value_HK']],
    'FC_vs_MC_HK': significant_FC_MC_HK[['Name', 'FC_to_MC_FoldChange_HK', 'FC_MC_p_value_HK']],
    'AF_vs_AM_NEG': significant_AF_AM_NEG[['Name', 'AF_to_AM_FoldChange_NEG', 'AF_AM_p_value_NEG']],
    'AF_vs_FC_NEG': significant_AF_FC_NEG[['Name', 'AF_to_FC_FoldChange_NEG', 'AF_FC_p_value_NEG']],
    'AM_vs_MC_NEG': significant_AM_MC_NEG[['Name', 'AM_to_MC_FoldChange_NEG', 'AM_MC_p_value_NEG']],
    'FC_vs_MC_NEG': significant_FC_MC_NEG[['Name', 'FC_to_MC_FoldChange_NEG', 'FC_MC_p_value_NEG']],
    'AF_vs_AM_POS': significant_AF_AM_POS[['Name', 'AF_to_AM_FoldChange_POS', 'AF_AM_p_value_POS']],
    'AF_vs_FC_POS': significant_AF_FC_POS[['Name', 'AF_to_FC_FoldChange_POS', 'AF_FC_p_value_POS']],
    'AM_vs_MC_POS': significant_AM_MC_POS[['Name', 'AM_to_MC_FoldChange_POS', 'AM_MC_p_value_POS']],
    'FC_vs_MC_POS': significant_FC_MC_POS[['Name', 'FC_to_MC_FoldChange_POS', 'FC_MC_p_value_POS']]
}

# Filter out empty DataFrames
significant_summary = {key: df for key, df in significant_summary.items() if not df.empty}

# Output the results summary
print("__________________________________________________________")
print("-------------Summary of Significant miRNAs----------------")
for key, df in significant_summary.items():
    print(f"\n{key}:")
    print(df)
print("==========================================================")


# Check the statistical distribution of p-values in the data
# Extract all p-value columns
p_value_columns = [col for col in grouped_data.columns if 'p_value' in col]
p_values = grouped_data[p_value_columns].values.flatten()

# Round p-values to the nearest 0.00
rounded_p_values = p_values.round(2)

# Calculate distribution
p_value_distribution = pd.Series(rounded_p_values).value_counts().sort_index().reset_index()

# Rename columns
p_value_distribution.columns = ['P-Value', 'Count']

# Output the distribution for validation
print("__________________________________________________________")
print("-----------P-Values Statistical Distribution--------------")
print(p_value_distribution)
print("==========================================================")

# File path for the p-values distribution results
p_value_distribution_file_path = './pvalues_miRNA_distribution.csv'

# Saving the data to a CSV file
p_value_distribution.to_csv(p_value_distribution_file_path, index=False)

# Concatenate all DataFrames into a single DataFrame
all_significant_miRNAs = pd.concat(significant_summary.values(), ignore_index=True)

# File path for the significant results
significant_results_summary_file_path = './significant_miRNA_results_summary.csv'

# Saving the data to a CSV file
all_significant_miRNAs.to_csv(significant_results_summary_file_path, index=False)

# File path for the grouped data with the fold changes and p-values
grouped_miRNA_folded_pvalues_file_path = './grouped_miRNA_folded_pvalues.csv'

# Saving the data to a CSV file
grouped_data.to_csv(grouped_miRNA_folded_pvalues_file_path, index=False)

exit()
