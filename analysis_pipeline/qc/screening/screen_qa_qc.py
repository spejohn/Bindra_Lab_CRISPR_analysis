import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn as sns

import os
import pandas as pd
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go
import logging

# Setup logger
logger = logging.getLogger(__name__)

def create_fill_QC_df(directory, output_dir, filename, overwrite = False, graph = False):
    final_output = Path(output_dir, filename)
    if final_output.exists() and overwrite is False:
        print(f'Output file already made, see {output_dir}\\{filename}. Change name or overwrite status.')
    else:
        
        # Initialize QC_df
        QC_list = {"First_Last": [], "Screen_ID":[], 
                    "neg_lfc_R2_MGK":[],
                     "neg_lfc_slope_MGK":[],
                     "neg_score_R2_MGK":[],
                     "neg_score_slope_MGK":[],
                     "pos_score_R2_MGK":[],
                     "pos_score_slope_MGK":[],
                     "neg_fdr_R2_MGK":[],
                     "neg_fdr_slope_MGK":[],
                     "-log(neg_fdr)_R2_MGK":[],
                     "-log(neg_fdr)_slope_MGK":[],
                     "pos_fdr_R2_MGK":[], 
                     "pos_fdr_slope_MGK":[],
                     "-log(pos_fdr)_R2_MGK":[],
                     "-log(pos_fdr)_slope_MGK":[],
                     "neg_rank_R2_MGK":[], 
                     "neg_rank_slope_MGK":[],
                     "num_R2_MGK":[],
                     "num_slope_MGK":[],

                     "normz_R2_DZ":[],
                     "normz_slope_DZ":[],
                     "fdr_synth_R2_DZ":[],
                     "fdr_synth_slope_DZ":[],
                     "-log(fdr_synth)_R2_DZ":[],
                     "-log(fdr_synth)_slope_DZ":[],
                     "fdr_supp_R2_DZ":[], 
                     "fdr_supp_slope_DZ":[],
                     "-log(fdr_supp)_R2_DZ":[],
                     "-log(fdr_supp)_slope_DZ":[],
                     "rank_synth_R2_DZ":[],
                     "rank_synth_slope_DZ":[],
                     "numobs_R2_DZ":[],
                     "numobs_slope_DZ":[]}
        
        if Path(output_dir + os.sep + "directory_df_output.csv").exists():
            directory_df = pd.read_csv(Path(output_dir + os.sep + "directory_df_output.csv"))
        else:
            directory_df = check_pairs(directory)
            #directory_df.to_csv(output_dir + os.sep + "directory_df_output.csv")

        # There were issues with overwriting/erasing indexes for "First_Last" and "Screen_ID" when filling QC_df/list with MGK/DZ_stats in the loop for some reason. I decided it would just be easiest to add all screen ID info at the beginning, index by both columns, then add MGK/DZ_stats values by index. I'm only using the provided column because the directory_df is "provided" and "generated" columns of paired _DZ/MGK.csv files so shouldn't be a difference b/w provided and generated anyways. P.S. the list has to be filled with empty strings because uneven arrays can't be converted to df later if values are empty. 
        for path in directory_df['provided']:
            screen_id = Path(path).parent.name
            first_last = Path(path).parent.parent.name
            new_rows = {"First_Last": [first_last], "Screen_ID": [screen_id]}
            for key in new_rows:
                QC_list[key].extend(new_rows[key])
            
            # Find the maximum length of lists in the dictionary
            max_length = max(len(lst) for lst in QC_list.values())

            # Add empty strings ('') to lists that are shorter than max_length
            for key in QC_list:
                QC_list[key] += [''] * (max_length - len(QC_list[key]))

        QC_df = pd.DataFrame(QC_list)
        QC_df.set_index(["First_Last", "Screen_ID"], inplace=True)

        for index, row in directory_df.iterrows():

            pDZ_df = None
            pMGK_df = None
            gDZ_df = None
            gMGK_df = None

            if '_pDZ' in row['provided']:
                screen_id = Path(row['provided']).parent.name
                first_last = Path(row['provided']).parent.parent.name

                DZ_output_dir = Path(os.path.join(output_dir, first_last, screen_id))

                pDZ_df = pd.read_csv(row['provided'])
                print(f"Reading {row['provided']} into a dataframe.")

            elif '_pMGK' in row['provided']:
                screen_id = Path(row['provided']).parent.name
                first_last = Path(row['provided']).parent.parent.name

                MGK_output_dir = Path(os.path.join(output_dir, first_last, screen_id))

                pMGK_df = pd.read_csv(row['provided'])
                print(f"Reading {row['provided']} into a dataframe.")

            else:
                print('Unrecognized file names. Please verify _pDZ/MGK is present.')

            if '_gDZ' in row['generated']:
                #gDZ_path = Path(row['generated'])

                gDZ_df = pd.read_csv(row['generated'])
                print(f"Reading {row['generated']} into a dataframe.")

            elif '_gMGK' in row['generated']:
                #gMGK_path = Path(row['generated'])

                gMGK_df = pd.read_csv(row['generated'])
                print(f"Reading {row['generated']} into a dataframe.")
            
            else:
                print('Unrecognized file names. Please verify _gDZ/MGK is present.')
            
            # Make all column names lowercase so case-insensitive.
            for i, df in enumerate([pDZ_df, pMGK_df, gDZ_df, gMGK_df]):
                if df is None:
                    continue
                else:
                    df.rename(columns=lambda x: x.lower(), inplace=True)

            # Output missing genes, calculate statistics, and graph paired files.
            if pDZ_df is not None and gDZ_df is not None:
                if len(pDZ_df['gene']) != len(gDZ_df['gene']):
                    # Find genes in gDZ_df that are not in pDZ_df
                    find_excluded_genes(pDZ_df, gDZ_df, DZ_output_dir)

                DZ_stats = None
                DZ_stats = calc_stats(pDZ_df, gDZ_df, 'DZ', DZ_output_dir, graph = graph)
                print(f"Calculated stats for pDZ and gDZ tables.")

                # Fill the statistics into the corresponding columns in QC_df
                if DZ_stats is not None:
                    for key, value in DZ_stats.items():
                        value = value.iloc[0] # Get the first value of the Series
                        if 'R2' in key:
                            QC_df.loc[(first_last, screen_id), f"{key}"] = value
                        elif 'slope' in key:
                            QC_df.loc[(first_last, screen_id), f"{key}"] = value


            if pMGK_df is not None and gMGK_df is not None:
                
                if len(pMGK_df['id']) != len(gMGK_df['id']):
                    # Find genes in gDZ_df that are not in pDZ_df
                    find_excluded_genes(pMGK_df, gMGK_df, MGK_output_dir)
                
                MGK_stats = None
                MGK_stats = calc_stats(pMGK_df, gMGK_df, 'MGK', MGK_output_dir, graph = graph)
                print(f"Calculated stats for pMGK and gMGK tables.")

                # Fill the statistics into the corresponding columns in QC_df
                if MGK_stats is not None:
                    for key, value in MGK_stats.items():
                        value = value.iloc[0] # Get the first value of the Series
                        if 'R2' in key:
                            QC_df.loc[(first_last, screen_id), f"{key}"] = value
                        elif 'slope' in key:
                            QC_df.loc[(first_last, screen_id), f"{key}"] = value
      
        output_path = Path(output_dir) / filename
        QC_df.to_csv(output_path)
        return QC_df


def check_pairs(directory: Path) -> pd.DataFrame:
    """
    Check for matching pairs of DrugZ and MAGeCK analysis files.
    
    Args:
        directory: Path to the directory containing analysis results
        
    Returns:
        DataFrame with columns for provided and generated file paths
    """
    dict_pairs = {'provided': [], 'generated': []}
    
    # Recursively iterate through provided directory
    for dir_path, subdirs, files in os.walk(directory):
        
        gDZ_files = []
        pDZ_files = []
        gMGK_files = []
        pMGK_files = []

        for file in files:
            file_path = os.path.join(dir_path, file)
            # Find files that match either p/gMGK or p/gDZ file patterns
            if file.endswith('_gDZ.csv'):
                gDZ_files.append(file_path)
                logger.debug(f'Found DrugZ generated file: {file}')
            elif file.endswith('_gMGK.csv'):
                gMGK_files.append(file_path)
                logger.debug(f'Found MAGeCK generated file: {file}')
            elif file.endswith('_pDZ.csv'):
                pDZ_files.append(file_path)
                logger.debug(f'Found DrugZ provided file: {file}')
            elif file.endswith('_pMGK.csv'):
                pMGK_files.append(file_path)
                logger.debug(f'Found MAGeCK provided file: {file}')
        
        # Match files by base name (without suffix)
        # For DrugZ files
        for gDZ_file in gDZ_files:
            g_basename = os.path.basename(gDZ_file).replace('_gDZ.csv', '')
            
            # Look for matching provided file
            for pDZ_file in pDZ_files:
                p_basename = os.path.basename(pDZ_file).replace('_pDZ.csv', '')
                
                if g_basename == p_basename:
                    dict_pairs['generated'].append(gDZ_file)
                    dict_pairs['provided'].append(pDZ_file)
                    logger.info(f'Matched DrugZ pair: {g_basename}')
                    break
        
        # For MAGeCK files
        for gMGK_file in gMGK_files:
            g_basename = os.path.basename(gMGK_file).replace('_gMGK.csv', '')
            
            # Look for matching provided file
            for pMGK_file in pMGK_files:
                p_basename = os.path.basename(pMGK_file).replace('_pMGK.csv', '')
                
                if g_basename == p_basename:
                    dict_pairs['generated'].append(gMGK_file)
                    dict_pairs['provided'].append(pMGK_file)
                    logger.info(f'Matched MAGeCK pair: {g_basename}')
                    break
    
    # Create DataFrame from dictionary
    df = pd.DataFrame(dict_pairs)
    
    # Log summary
    if not df.empty:
        logger.info(f'Found {len(df)} matching file pairs')
    else:
        logger.warning('No matching file pairs found')
        
    return df

def find_excluded_genes(p_df: pd.DataFrame, g_df: pd.DataFrame, output_dir: Path):
    #p_df should be 'provided' dataframe (pDZ/MGK) and g_df should be 'generated' dataframe (gDZ/MGK)
    if 'gene' in p_df.columns and 'gene' in g_df.columns:
        gene_column = 'gene'
    elif 'id' in p_df.columns and 'id' in g_df.columns:
        gene_column = 'id'
    else:
        raise ValueError("Both dataframes must contain a common name for the gene column.")

    # Find genes in generated_df that are not in provded_df
    excluded_genes = g_df[~g_df[gene_column].isin(p_df[gene_column])]

    # Save the result to a new CSV file
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        excluded_genes.to_csv(output_dir / '_genes_from_gfile_not_in_pfile.csv', index=False)
        print(f"Excluded genes have been saved to {output_dir}")

def calc_stats(p_df: pd.DataFrame, g_df: pd.DataFrame, analysis_type: str, output_dir: Path, graph = False):
    if analysis_type == 'DZ':
        join_column = 'gene'
    elif analysis_type == 'MGK':
        join_column = 'id'
    else:
        raise ValueError("Invalid match_column. Please specify 'DZ' or 'MGK'.")
    
    p_df.columns = p_df.columns.str.replace('|', '_')
    g_df.columns = g_df.columns.str.replace('|', '_')

    # Match rows by the determined column
    p_df = p_df.set_index(join_column)
    g_df = g_df.set_index(join_column)
    p_df, g_df = p_df.align(g_df, join='inner')
    p_df = p_df.reset_index()
    g_df = g_df.reset_index()

    p_df = apply_neg_log(p_df)
    g_df = apply_neg_log(g_df)

    data_columns = ['normz', 'fdr_synth', 'fdr_supp', 'rank_synth', 'neg_lfc', 'neg_score', 'pos_score', 'neg_fdr', 'pos_fdr', 'neg_rank', '-log(fdr_synth)', '-log(fdr_supp)', '-log(neg_fdr)', '-log(pos_fdr)', 'numobs', 'num']
    
    R_list = {}
    slope_list = {}
    
    for column in data_columns:
        if column in p_df.columns and column in g_df.columns:
            X = p_df[column].values.reshape(-1, 1)
            Y = g_df[column].values
            reg = LinearRegression().fit(X, Y)
            r_squared = round(reg.score(X, Y), 3)
            slope = round(reg.coef_[0], 3)

            R_list[column + (f"_R2_{analysis_type}")] = r_squared # Replace the list with a new list containing one value
            slope_list[column + (f"_slope_{analysis_type}")] = slope # Replace the list with a new list containing one value

            R_slope_list = pd.DataFrame({**R_list, **slope_list}, index = [0])

            if graph == True:

                # Create a scatter plot
                fig = px.scatter(p_df, x = X.squeeze(), y = Y, title = f'{column}: R-squared = {r_squared}, Slope = {slope}', hover_name = join_column)
                fig.update_layout(xaxis_title = "provided values", yaxis_title = "generated values") 
            
                # Add the regression line
                fig.add_trace(go.Scatter(x = X.squeeze(), y = reg.predict(X), mode='lines', name = 'Regression Line'))

                if not output_dir.exists():
                    output_dir.mkdir(parents=True)

                fig.write_html(
                Path(output_dir)
                / f"{column}_plot.html")

                print(f"Linear regression plot for {column} saved in {output_dir}.")

    return R_slope_list

def apply_neg_log(df):
    # Select columns whose names contain 'fdr'
    fdr_columns = [col for col in df.columns if 'fdr' in col]

    # Apply the negative logarithm to the selected columns and append them to the DataFrame
    for col in fdr_columns:
        df[f"-log({col})"] = df[col].apply(lambda x: -np.log(x))

    return df



data_path = "C:\\Users\\spejo\\Box\\DDR_Brow\\Datasets\\4_Completed_and_uploaded"
output_folder = "C:\\Users\\spejo\\Box\\DDR_Brow\\QA_QC"
test_folder = "C:\\Users\\spejo\\Documents\\BI_Pipeline\\output"

create_fill_QC_df(data_path, output_folder, filename = "7-19-24_screen_QA_analysis.csv", graph = False, overwrite = False)