# miRNA_data_analysis
 Analyze and Visualize Nanostring miRNA sample data


# miRNA_analyzer.py
This code is a comprehensive script for performing a variety of tasks with your NanoString data, including:

- Loading the data and performing initial checks.
- Identifying housekeeping, negative, and positive control miRNAs from a list of externally calculated stable miRNAs.
- Normalizing expression levels to housekeeping, negative, and positive controls.
- Calculating average expression and fold changes between groups.
- Performing t-tests to assess statistical differences between groups.
- Creates a summary of "significant" results as defined by the fold change and p-value threshhold variables (defaults: foldchange >= 2, p-value <= 0.1)
- Saving the final processed data to CSV files, one for the grouped data and one for the significant results summary data.

Before you run this script, please ensure that:

1. Your dataset is correctly formatted and includes 'Code Class', 'Accession', 'Name', and the expression columns with appropriate headers.
2. The expression columns for the samples should be prefixed appropriately so that the script can differentiate between different groups (AF, AM, FC, MC).
3. You have installed the required Python packages (`pandas` and `scipy`) in your environment.

## Script Overview
**Filename**: miRNA_analyzer.py  
**Description**: This script analyzes microRNA (miRNA) expression data to identify significant changes in expression levels between different sample groups. It performs normalization using housekeeping, negative, and positive miRNAs, calculates fold changes, performs t-tests, and identifies significant miRNAs based on fold change and p-value thresholds.

## Dependencies
- pandas: For data manipulation and analysis.
- scipy.stats: Specifically `ttest_ind` for performing independent t-tests.

## Script Workflow
1. **Data Loading and Initial Processing**:
   - Loads miRNA data from the provided CSV file into a pandas DataFrame.
   - Sets predefined lists of stable reference miRNAs, including housekeeping, negative, and positive miRNAs.

2. **Normalization**:
   - It identifies housekeeping, negative, and positive miRNAs in the input data set.
   - Normalizes the expression data using average expressions of these miRNAs.

3. **Grouping and Fold Change Calculation**:
   - The script groups the data by sample types (AF, AM, FC, MC) for all normalization methods.
   - Calculates the mean expression values for each group.
   - Computes fold changes between different sample groups.

4. **Statistical Testing**:
   - Performs t-tests to compare miRNA expressions between different sample groups.
   - Calculates p-values for these comparisons.

5. **Identification of Significant miRNAs**:
   - Selects miRNAs with significant fold changes (threshold > 2) and low p-values (threshold < 0.1).
   - Summarizes these significant miRNAs in a dictionary.

6. **Results Aggregation and Export**:
   - Aggregates significant miRNAs from all comparisons into a single DataFrame.
   - Exports this aggregated data and the full grouped dataset with fold changes and p-values to CSV files.

## Functions
### `perform_ttest(data_normalized, group1, group2)`
Calculates p-values using t-tests for normalized groups.

**Parameters**:
- `data_normalized`: DataFrame containing grouped & normalized data (e.g. housekeeping_normalized_miRNAs)
- `group1`, `group2`: Column names of the groups to be compared.

**Returns**: List of p-values for each miRNA.

### `select_significant_miRNAs(data_grouped, fold_change_column, p_value_column)`
Selects miRNAs that meet specified fold change and p-value thresholds.

**Parameters**:
- `data_grouped`: DataFrame containing miRNA expression data.
- `fold_change_column`: Column name for fold change values.
- `p_value_column`: Column name for p-values.

**Returns**: DataFrame with significant miRNAs based on the specified criteria.

## Usage
Run the script in a Python environment where all dependencies are installed. Ensure the input CSV file path is correctly specified.

