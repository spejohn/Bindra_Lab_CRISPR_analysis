import os
import re
import subprocess
from pathlib import Path, PurePosixPath
from typing import Optional, Union
import logging
from datetime import datetime

import docker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
from plotly import graph_objects as go
from scipy.stats import pearsonr
from sklearn.metrics import auc, precision_recall_curve, roc_curve
from profilehooks import profile

DOCKER_CLIENT = docker.from_env()

# Configure logging
timestamp = datetime.now().strftime('%m-%d-%y_%H-%M')  # MM-DD-YY_HH-MM format
log_file = rf'C:\Users\spejo\Documents\2_CRISPR_analysis_test_output\validation_errors_{timestamp}.log'

logging.basicConfig(
    filename=log_file, 
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def generate_sample_sheet(fastq_dir: str) -> pd.DataFrame:
    """
    Generate a sample sheet DataFrame from FASTQ files in a directory.
    """
    fastq_files = []
    fastq_dir_path = Path(fastq_dir)
    
    # Find all FASTQ files recursively
    for file in fastq_dir_path.rglob('*'):
        if file.suffix.lower() in ['.fastq', '.fq']:
            fastq_files.append(file)
    
    if not fastq_files:
        raise FileNotFoundError(f"No FASTQ files found in {fastq_dir}")
    
    # Create sample sheet with absolute paths
    sample_sheet = pd.DataFrame({
        'sample_name': [f.stem for f in fastq_files],
        'fastq_path': [str(f.absolute()) for f in fastq_files]
    })
    
    logging.info(f"Found {len(sample_sheet)} FASTQ files in {fastq_dir}")
    logging.debug(f"Sample sheet contents:\n{sample_sheet.to_string()}")
    
    return sample_sheet


def mageck_count(sample_sheet: pd.DataFrame, library_file: str, output_dir: str):
    """Run MAGeCK count to generate read counts from FASTQ files."""
    image = "samburger/mageck"
    
    # Convert all paths to absolute Path objects
    library_path = Path(library_file).absolute()
    output_path = Path(output_dir).absolute()
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Get experiment name from FASTQ directory
    fastq_dir = Path(os.path.dirname(sample_sheet['fastq_path'].iloc[0]))
    experiment_name = fastq_dir.name
    
    # Initialize volumes dictionary with just library and output paths
    volumes = {
        str(Path(library_file).parent): {'bind': '/library', 'mode': 'ro'},
        str(output_path): {'bind': '/output', 'mode': 'rw'},
    }
    
    # Add each unique FASTQ directory to volumes
    unique_fastq_dirs = {str(Path(path).parent) for path in sample_sheet['fastq_path']}
    for idx, fastq_dir in enumerate(unique_fastq_dirs):
        volumes[fastq_dir] = {'bind': f'/fastq/dir_{idx}', 'mode': 'ro'}
    
    try:
        # Build mageck count command
        count_command = [
            "mageck",
            "count",
            "-l", f"/library/{Path(library_file).name}",
            "-n", f"/output/{experiment_name}_RC",
            "--trim-5", "0",
            "--sample-label", ",".join(sample_sheet["sample_name"])
        ]

        # Build list of FASTQ paths in container
        fastq_paths = []
        for idx, row in sample_sheet.iterrows():
            fastq_dir_idx = list(unique_fastq_dirs).index(str(Path(row['fastq_path']).parent))
            container_path = f'/fastq/dir_{fastq_dir_idx}/{Path(row["fastq_path"]).name}'
            fastq_paths.append(container_path)
        
        # Add --fastq flag and all FASTQ files
        count_command.extend(["--fastq"] + fastq_paths)

        command_str = " ".join(map(str, count_command))
        logging.info(f"Running MAGeCK count command: {command_str}")
        logging.info(f"Volume mappings:")
        for host_path, container_info in volumes.items():
            logging.info(f"  {host_path} -> {container_info['bind']} ({container_info['mode']})")
        
        # Run container with explicit working directory
        container = DOCKER_CLIENT.containers.run(
            image,
            f'"{command_str}"',
            volumes=volumes,
            remove=True,
            detach=True,
            stdout=True,
            stderr=True
        )
        
        # Stream logs
        for line in container.logs(stream=True, follow=True):
            log_text = line.decode().strip()
            print(log_text)  # Print to console for immediate feedback
            logging.info(f"Container output: {log_text}")
        
        # Check container exit status
        result = container.wait()
        if result['StatusCode'] != 0:
            error_msg = f"Container exit code: {result['StatusCode']}\n"
            error_msg += f"Error message: {result.get('Error', 'No error message available')}"
            logging.error(error_msg)
            raise RuntimeError(f"MAGeCK count failed: {error_msg}")
        
        # List output directory contents for debugging
        logging.info(f"Output directory contents: {list(output_path.glob('*'))}")
        
        # Verify output files exist and have content
        count_file = output_path / f"{experiment_name}_RC.count.txt"
        if not count_file.exists():
            raise FileNotFoundError(f"Expected count file not found: {count_file}")
        
        # Read and verify count table
        try:
            count_data = pd.read_csv(count_file, sep='\t')
            if count_data.empty:
                raise ValueError(f"Count table is empty: {count_file}")
            logging.info(f"Count table has {len(count_data)} rows and {len(count_data.columns)} columns")
        except Exception as e:
            logging.error(f"Error reading count table: {str(e)}")
            raise
        
        # Create CSV version with absolute path
        new_output_file = output_path / f"{experiment_name}_RC.csv"
        count_data.to_csv(new_output_file, index=False)
        logging.info(f"Generated CSV count table: {new_output_file}")
        
        # Add this after running the container:
        container_logs = container.logs().decode()
        if "error" in container_logs.lower() or "exception" in container_logs.lower():
            logging.error(f"Container logs indicate error:\n{container_logs}")
        
        # Check container filesystem
        try:
            output_ls = container.exec_run("ls -la /output").output.decode()
            logging.info(f"Container output directory contents:\n{output_ls}")
        except Exception as e:
            logging.error(f"Could not check container filesystem: {str(e)}")
        
        return str(new_output_file)

    except Exception as e:
        logging.error(f"Error running MAGeCK count: {str(e)}")
        logging.error(f"Command: {command_str}")
        logging.error(f"Volumes: {volumes}")
        logging.error(f"Sample sheet contents:\n{sample_sheet.to_string()}")
        raise


def mageck_test(contrasts: str, input_file: str, output_dir: str):
    """
    Run MAGeCK test analysis.
    Returns a list of output file paths.
    """
    image = "samburger/mageck"
    output_files = []
    
    # Convert paths to absolute Path objects
    input_path = Path(input_file).absolute()
    output_path = Path(output_dir).absolute()
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Match the working volume format exactly
    volumes = [
        f"{str(Path(os.getcwd()).absolute())}:/work:rw",
        f"{str(output_path)}:/output:rw"
    ]
    
    logging.info(f"Input file: {input_path}")
    logging.info(f"Output directory: {output_path}")
    logging.info(f"Volume mappings: {volumes}")
    
    # Read contrasts file
    contrasts_df = pd.read_csv(Path(contrasts), sep='\t')
    
    # Process each contrast
    for line in contrasts_df.itertuples():
        try:
            # Handle input file based on extension
            if input_path.suffix.lower() == '.csv':
                tsv_name = f"{input_path.stem.replace('_RC', '')}.counts.txt"
                host_tsv = output_path / tsv_name
                count_data = pd.read_csv(input_path)
                count_data.to_csv(host_tsv, sep='\t', index=False)
                logging.info(f"Converted {input_path} to tab-separated format: {host_tsv}")
                container_input = f"/output/{tsv_name}"
            else:
                container_input = f"/work/{input_path.name}"
            
            # Build command
            test_command = [
                "mageck",
                "test",
                "--count-table",
                container_input,
                "--norm-method",
                "median",
                "--adjust-method",
                "fdr",
                "--output-prefix",
                f"/output/{line.contrast}",
                "--keep-tmp",
                "-t", line.treatment,
                "-c", line.control
            ]
            
            command_str = " ".join(map(str, test_command))
            logging.info(f"Running MAGeCK command for contrast {line.contrast}:")
            logging.info(f"Command: {command_str}")
            
            # Match the working container run format exactly
            shell_cmd = f"docker run -it"
            for vol in volumes:
                shell_cmd += f' -v "{vol}"'
            shell_cmd += f' samburger/mageck "{command_str}"'
            logging.info(f"Shell command: {shell_cmd}")
            
            container = DOCKER_CLIENT.containers.run(
                image,
                f'"{command_str}"',  # Match the working format with quotes
                volumes=volumes,
                remove=False,
                detach=True
            )

            # Stream logs
            for log_line in container.logs(stream=True, follow=True):
                log_text = log_line.decode().strip()
                print(log_text)
                logging.info(f"Container output: {log_text}")
            
            result = container.wait()
            if result['StatusCode'] != 0:
                raise RuntimeError(f"MAGeCK test failed with status code {result['StatusCode']}")
                
            # Verify output files exist
            expected_files = [
                output_path / f"{line.contrast}.gene_summary.txt",
                output_path / f"{line.contrast}.sgrna_summary.txt"
            ]
            
            for file in expected_files:
                logging.info(f"Checking for output file: {file}")
                if not file.exists():
                    logging.error(f"File not found: {file}")
                    logging.error(f"Directory contents: {list(output_path.glob('*'))}")
                    raise FileNotFoundError(f"Expected output file not found: {file}")
                
                # Verify file has content
                try:
                    df = pd.read_csv(file, sep='\t')
                    if df.empty:
                        raise ValueError(f"Output file is empty: {file}")
                    logging.info(f"Generated {file.name} with {len(df)} rows")
                    output_files.append(str(file))
                except Exception as e:
                    logging.error(f"Error reading output file {file}: {str(e)}")
                    raise
            
            # Convert gene summary to CSV
            gene_summary = expected_files[0]
            csv_file = output_path / f"{line.contrast}_gMGK.csv"
            pd.read_csv(gene_summary, sep="\t").to_csv(csv_file, index=False)
            output_files.append(str(csv_file))
            logging.info(f"Generated CSV output: {csv_file}")
            
        except Exception as e:
            logging.error(f"Error running MAGeCK test for contrast {line.contrast}")
            logging.error(f"Command: {command_str}")
            logging.error(f"Volumes: {volumes}")
            logging.error(f"Input file: {input_file}")
            logging.error(f"Error: {str(e)}")
            raise
                
    return output_files


def run_drugz(contrasts: str, input_file: str, output_dir: str):
    """
    Run DrugZ analysis.
    Returns a list of output file paths.
    """
    drugz_path = Path(os.getcwd()) / "drugz" / "drugz.py"
    if not drugz_path.exists():
        raise ModuleNotFoundError(
            "Drugz not found. Ensure you have run `git submodule update --init --recursive`."
        )

    # Convert paths to absolute Path objects
    input_path = Path(input_file).absolute()
    output_path = Path(output_dir).absolute()
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Read contrasts file
    contrasts_df = pd.read_csv(Path(contrasts), sep='\t')
    output_files = []
    
    # Process each contrast
    for c in contrasts_df.itertuples():
        try:
            # Handle input file based on extension
            counts_file = output_path / f"{c.contrast}_counts.txt"
            
            # Read and verify input file
            try:
                if input_path.suffix.lower() == '.csv':
                    count_data = pd.read_csv(input_path)
                else:
                    count_data = pd.read_csv(input_path, sep='\t')
                
                if count_data.empty:
                    raise ValueError(f"Input count table is empty: {input_path}")
                    
                # Always write as TSV for DrugZ
                count_data.to_csv(counts_file, sep='\t', index=False)
                logging.info(f"Prepared count data for DrugZ at {counts_file} with {len(count_data)} rows")
            except Exception as e:
                logging.error(f"Error processing input file: {str(e)}")
                raise
            
            output_file = output_path / f"{c.contrast}_DZ_output.tsv"
            drugz_command = [
                "python",
                str(drugz_path.absolute()),
                "-i",
                str(counts_file.absolute()),
                "-o",
                str(output_file.absolute()),
                "-c",
                c.control,
                "-x",
                c.treatment,
                "-unpaired"
            ]
            
            logging.info(f"Running DrugZ command for contrast {c.contrast}: {' '.join(drugz_command)}")
            
            # Capture both stdout and stderr
            result = subprocess.run(
                drugz_command, 
                check=True,
                capture_output=True,
                text=True
            )
            
            # Log the output
            if result.stdout:
                logging.info(f"DrugZ stdout:\n{result.stdout}")
            if result.stderr:
                logging.warning(f"DrugZ stderr:\n{result.stderr}")
            
            logging.info(f"Successfully ran DrugZ command for contrast {c.contrast}")
            
            # Verify output file exists and has content
            if not output_file.exists():
                raise FileNotFoundError(f"Expected DrugZ output file not found: {output_file}")
            
            try:
                drugz_data = pd.read_csv(output_file, sep="\t")
                if drugz_data.empty:
                    raise ValueError(f"DrugZ output file is empty: {output_file}")
                logging.info(f"Generated DrugZ output with {len(drugz_data)} rows")
                output_files.append(str(output_file))
            except Exception as e:
                logging.error(f"Error reading DrugZ output: {str(e)}")
                raise

            # Create new .csv for readability/downstream QA_QC
            csv_file = output_path / f"{c.contrast}_gDZ.csv"
            drugz_data.to_csv(csv_file, index=False)
            output_files.append(str(csv_file))
            logging.info(f"Generated CSV output: {csv_file}")
            
        except subprocess.CalledProcessError as e:
            logging.error(f"DrugZ failed for contrast {c.contrast}")
            logging.error(f"Command: {' '.join(drugz_command)}")
            logging.error(f"Error output: {e.stderr if e.stderr else 'No error output'}")
            raise
        except Exception as e:
            logging.error(f"Error processing DrugZ output for contrast {c.contrast}: {str(e)}")
            raise

    return output_files


def plot_QA(
    screen_title: str,
    ess_genes: dict[str, list[str]],
    noness_genes: Optional[dict[str, list[str]]] = None,
):
    # TODO: add guide- and gene-level count correlation between replicates
    # TODO: make figures dir if does not exist
    # Load mageck sgRNA counts
    sgrna_counts = pd.read_csv(
        Path("data") / "mageck" / screen_title / f"{screen_title}.count.txt", sep="\t"
    )
    plot_QA_sgRNA_corr(screen_title, sgrna_counts)

    # Load mageck gene LFCs
    dfs = []
    for filename in Path(f"data/mageck/{screen_title}").iterdir():
        if "gene_summary" in filename.name:
            contrast = filename.name[
                filename.name.find("_") + 1 : filename.name.find(".gene_summary")
            ]
            df = pd.read_csv(filename, sep="\t")
            df["contrast"] = contrast
            dfs.append(df)
    mageck_results = pd.concat(dfs, ignore_index=True)
    # drop nontargeting guides
    mageck_results = mageck_results[mageck_results["id"] != "Non_Targeting_Human_CRko"]
    for genelist in ess_genes:
        if genelist not in noness_genes:
            # Label all genes not in essential list as non-essential.
            mageck_results["essential"] = mageck_results["id"].apply(
                lambda x: 1 if x in ess_genes[genelist] else 0
            )
            plot_QA_ROC(screen_title, mageck_results, genelist)
            plot_QA_PRC(screen_title, mageck_results, genelist)
        else:
            # Limit to only genes in essential/nonessential lists
            mageck_results = mageck_results[
                (mageck_results["id"].isin(ess_genes[genelist]))
                | (mageck_results["id"].isin(noness_genes[genelist]))
            ]
            mageck_results["essential"] = mageck_results["id"].apply(
                lambda x: 1 if x in ess_genes[genelist] else 0
            )
            plot_QA_ROC(screen_title, mageck_results, genelist)
            plot_QA_PRC(screen_title, mageck_results, genelist)


def plot_QA_sgRNA_corr(screen_title: str, counts: pd.DataFrame):
    print("Generating sgRNA replicate correlation plots...")
    # Create a dictionary to hold the replicates
    replicate_dict = {}
    # Iterate over the columns and populate the dictionary with replicates
    for column in counts.columns:
        if "Rep" in column:
            # Extract the sample name and replicate number
            match = re.findall("Sample[_-](.*)[_-]Rep", column)
            if match:
                sample_name = match[0]
                replicate_number = re.findall(r"Rep[_-]?(\d+)", column)[0]
                # Add the replicate to the dictionary
                replicate_dict.setdefault(sample_name, []).append(
                    (replicate_number, column)
                )
    # Sort the replicates and pair them up
    replicates = []
    for sample, reps in replicate_dict.items():
        # Sort by replicate number
        sorted_reps = sorted(reps, key=lambda x: int(x[0]))
        # Pair up the replicates, assuming they are in order
        for i in range(0, len(sorted_reps) - 1, 2):
            replicates.append((sorted_reps[i][1], sorted_reps[i + 1][1]))
    for rep1, rep2 in replicates:
        plt.figure()
        g = sns.jointplot(
            x=rep1,
            y=rep2,
            data=counts,
            kind="reg",
            joint_kws={"line_kws": {"color": "black", "linewidth": 1}},
        )

        rep_guide_pearsonr = pearsonr(counts[rep1], counts[rep2]).statistic
        sample_title = re.findall("Sample[_-](.*)[_-]Rep", rep1)[0]
        plt.suptitle(
            f"{sample_title}: " + r"$\rho=$" + f"{rep_guide_pearsonr:.3f}", y=1.05
        )
        plt.xlabel("Rep 1 guide counts")
        plt.ylabel("Rep 2 guide counts")
        plt.savefig(
            Path("figures")
            / screen_title
            / f"QA_sgRNA_count_corr_{screen_title}_{sample_title}.png"
        )
        print(f'"{screen_title}","{sample_title}",{rep_guide_pearsonr},"pearson_r"')


def plot_QA_ROC(
    screen_title: str,
    mageck_results: pd.DataFrame,
    genelist: str,
):
    print("Generating ROC plots...")
    plt.figure()
    for contrast in mageck_results["contrast"].unique():
        subset = mageck_results[mageck_results["contrast"] == contrast].sort_values(
            "neg|lfc", ascending=False
        )
        fpr, tpr, _ = roc_curve(subset["essential"], -subset["neg|lfc"])
        roc_auc = auc(fpr, tpr)
        print(f'"{screen_title}","{contrast}",{roc_auc},"ROC_AUC"')
        plt.plot(fpr, tpr, lw=2, label=f"{contrast} (AUC = {roc_auc:0.2f})")
    plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC - {screen_title}")
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.savefig(
        Path("figures") / screen_title / f"QA_ROC_{screen_title}_{genelist}.png",
        bbox_inches="tight",
    )
    # Make plotly chart for debug with separate lines for each contrast
    fig = go.Figure()
    for contrast in mageck_results["contrast"].unique():
        subset = mageck_results[mageck_results["contrast"] == contrast].sort_values(
            "neg|lfc", ascending=False
        )
        fpr, tpr, threshold = roc_curve(subset["essential"], -subset["neg|lfc"])
        roc_auc = auc(fpr, tpr)
        roc_data = pd.DataFrame(
            {
                "threshold": threshold,
                "fpr": fpr,
                "tpr": tpr,
            }
        )
        # count ess genes at each threshold
        roc_data["n_genes_threshold"] = roc_data["threshold"].apply(
            lambda x: (-subset["neg|lfc"] > x).sum()
        )
        fig.add_trace(
            go.Scatter(
                x=roc_data.fpr,
                y=roc_data.tpr,
                mode="lines",
                name=f"{contrast} (AUC = {roc_auc:0.2f})",
                customdata=roc_data[["threshold", "n_genes_threshold"]],
                hovertemplate="<br>".join(
                    [
                        "LFC Threshold: %{customdata[0]}",
                        "TPR: %{y}",
                        "FPR: %{x}",
                        "# Genes > Threshold: %{customdata[1]}",
                    ]
                ),
            )
        )

    fig.update_layout(
        title=f"ROC Curve - {screen_title}",
        xaxis_title="False Positive Rate",
        yaxis_title="True Positive Rate",
        legend_title="Contrasts",
    )
    fig.write_html(
        Path("figures") / screen_title / f"QA_ROC_{screen_title}_{genelist}.html"
    )


def plot_QA_PRC(screen_title: str, mageck_results: pd.DataFrame, genelist: str):
    print("Generating Precision-Recall plots...")
    plt.figure()
    for contrast in mageck_results["contrast"].unique():
        subset = mageck_results[mageck_results["contrast"] == contrast].sort_values(
            "neg|lfc", ascending=False
        )
        precision, recall, _ = precision_recall_curve(
            subset["essential"], -subset["neg|lfc"]
        )
        pr_auc = auc(recall, precision)
        print(f'"{screen_title}","{contrast}",{pr_auc},"PR_AUC"')
        plt.plot(recall, precision, lw=2, label=f"{contrast} (AUC = {pr_auc:0.2f})")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"P/R - {screen_title}")
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.savefig(
        Path("figures") / screen_title / f"QA_PRC_{screen_title}_{genelist}.png",
        bbox_inches="tight",
    )
    # Make plotly chart for debug with separate lines for each contrast
    fig = go.Figure()
    for contrast in mageck_results["contrast"].unique():
        subset = mageck_results[mageck_results["contrast"] == contrast].sort_values(
            "neg|lfc", ascending=False
        )
        precision, recall, threshold = precision_recall_curve(
            subset["essential"], -subset["neg|lfc"]
        )
        pr_auc = auc(recall, precision)
        pr_data = pd.DataFrame(
            {
                "threshold": threshold,
                "precision": precision[:-1],
                "recall": recall[:-1],
            }
        )
        # count ess genes at each threshold
        pr_data["n_genes_threshold"] = pr_data["threshold"].apply(
            lambda x: (-subset["neg|lfc"] > x).sum()
        )
        fig.add_trace(
            go.Scatter(
                x=pr_data.recall,
                y=pr_data.precision,
                mode="lines",
                name=f"{contrast} (AUC = {pr_auc:0.2f})",
                customdata=pr_data[["threshold", "n_genes_threshold"]],
                hovertemplate="<br>".join(
                    [
                        "LFC Threshold: %{customdata[0]}",
                        "Precision: %{y}",
                        "Recall: %{x}",
                        "# Genes > Threshold: %{customdata[1]}",
                    ]
                ),
            )
        )

    fig.update_layout(
        title=f"Precision-Recall Curve - {screen_title}",
        xaxis_title="Recall",
        yaxis_title="Precision",
        legend_title="Contrasts",
    )
    fig.write_html(
        Path("figures") / screen_title / f"QA_PRC_{screen_title}_{genelist}.html"
    )


def plot_hits_mageck():
    # TODO: volcano plots and essentiality rank plots for mageck outputs
    pass


def plot_hits_drugz(screen_title: str, contrasts: pd.DataFrame, output_dir: str):
    # Load each drugz output file and combine to a single dataframe
    dfs = []
    for c in contrasts.itertuples():
        df = pd.read_csv(
            Path(output_dir)
            / "drugz"
            / screen_title
            / f"{screen_title}_{c.contrast}.tsv",
            sep="\t",
        )
        df["Contrast"] = c.contrast
        dfs.append(df)
    drugz_results = pd.concat(dfs, ignore_index=True)
    # drop nontargeting
    drugz_results = drugz_results[drugz_results["GENE"] != "Non_Targeting_Human_CRko"]
    drugz_results = drugz_results[
        ["Contrast", "GENE"] + list(drugz_results.columns[1:-1])
    ]
    drugz_results = drugz_results.sort_values(
        by=["Contrast", "normZ"], ascending=[True, True]
    )
    drugz_results = drugz_results.reset_index(drop=True)
    drugz_results.to_csv(
        Path(output_dir) / "drugz" / screen_title / f"drugz_results_{screen_title}.csv",
        index=False,
    )
    # Generate figures
    for contrast in drugz_results["Contrast"].unique():
        df = drugz_results.loc[drugz_results["Contrast"] == contrast]
        fig = px.scatter(
            df,
            x="rank_synth",
            y="normZ",
            hover_name="GENE",
            labels={"rank_synth": "Gene Rank"},
            title=contrast,
        )
        fig.write_html(
            Path(output_dir)
            / "drugz"
            / screen_title
            / f"{screen_title}_{contrast}_drugz_plot.html"
        )
    return

def validate_input_tables(
    contrasts: str,
    sample_sheet: Optional[pd.DataFrame] = None,
    counts: Optional[str] = None,
) -> bool:
    
    contrasts_df = pd.read_csv(Path(contrasts), sep='\t')
    try:
        assert set(contrasts_df.columns) == {"contrast", "treatment", "control"}
    except AssertionError:
        error_message = f"Error in {contrasts} column names - verify 'contrast', 'treatment', 'control' in table"
        print(error_message)
        logging.error(error_message)
        return False
    
    if counts is not None:
        try:
            # First try comma separator
            try:
                counts_df = pd.read_csv(Path(counts))
                logging.info(f"Successfully read {counts} with comma separator")
            except:
                # If that fails, try tab separator
                try:
                    counts_df = pd.read_csv(Path(counts), sep='\t')
                    logging.info(f"Successfully read {counts} with tab separator")
                except Exception as e:
                    logging.error(f"Failed to read {counts} with both comma and tab separators: {str(e)}")
                    return False
            
            # Log original columns
            logging.info(f"Original columns in {counts}: {counts_df.columns.tolist()}")
            
            # Convert column names to lowercase and strip whitespace
            columns_lower = {col.lower().strip() for col in counts_df.columns}
            required_columns = {'sgrna', 'gene'}
            
            # Log processed columns
            logging.info(f"Processed lowercase columns: {sorted(list(columns_lower))}")
            
            # Check for similar column names that might be variants
            for col in counts_df.columns:
                if 'guide' in col.lower() or 'sgrna' in col.lower():
                    logging.info(f"Found potential guide column: {col}")
                if 'gene' in col.lower():
                    logging.info(f"Found potential gene column: {col}")
            
            if not required_columns.issubset(columns_lower):
                missing_cols = required_columns - columns_lower
                error_message = (
                    f"{counts} is missing required columns: {missing_cols}.\n"
                    f"Found columns (original): {counts_df.columns.tolist()}\n"
                    f"Found columns (lowercase): {sorted(list(columns_lower))}\n"
                    f"First few rows:\n{counts_df.head().to_string()}"
                )
                print(error_message)
                logging.error(error_message)
                return False
                
            logging.info(f"Successfully validated count table {counts}")
            
        except Exception as e:
            error_message = f"Error reading count table {counts}: {str(e)}"
            print(error_message)
            logging.error(error_message)
            return False

    if sample_sheet is not None:
        try:
            assert set(sample_sheet.columns) == {"sample_name", "fastq_path"}
        except AssertionError:
            error_message = "Sample sheet columns are expected to be exactly: 'sample_name' and 'fastq_path'"
            print(error_message)
            logging.error(error_message)
            return False
    return True

def reverse_dir(input_dir:str, root:str, output_dir:str):
    # Use regex to find and remove the input_dir path from root
    if input_dir in root:
        dir_str = re.sub(f'^{re.escape(input_dir)}', '', root)
    else:
        print(f"Input directory string not found for {root} file.")

    # Combine the paths to mirror the input directory
    out_path = Path(output_dir) / Path(dir_str.lstrip('\\'))

    # Create dir if doesn't exist
    if not out_path.exists:
        os.makedirs(out_path, parents=True)

    return out_path

def check_existing_files(input_dir: str, output_dir: str):
    input_path = Path(input_dir)
    analyzed_files = set()  # Track specific files that have been analyzed

    # Iterate through input_dir for files ending in "_RC" or ".fastq"
    for root, dirs, files in os.walk(input_path):
        if dirs:
            continue
        for file in files:
            # Skip xlsx files
            if file.endswith('.xlsx'):
                continue
                
            if "_rc" in file.lower() or file.endswith(('.fq', '.fastq')):
                existing_output = reverse_dir(input_dir=input_dir, root=root, output_dir=output_dir)
                
                MGK_output_exists = any(existing_output.glob("*_gMGK.csv"))
                DZ_output_exists = any(existing_output.glob("*_gDZ.csv"))

                # Check if the parent directory exists in the output directory
                if MGK_output_exists or DZ_output_exists:
                    print(f"\nAnalysis for {file} exists in {existing_output}")
                    # Add the full path of the input file to analyzed_files
                    analyzed_files.add(str(Path(root) / file))
                    
            elif '_cnttbl' in file.lower():
                continue
            else:
                logging.debug(f"Skipping file {file} - not a recognized input file type")
    
    if not analyzed_files:
        logging.info(f"No previously analyzed files found in {input_dir}")
        
    return analyzed_files  # Return set of specific file paths that have been analyzed

def assign_file_path(root:str, file_str:str, fastq_dirs: set = None):
    if fastq_dirs is None:
        fastq_dirs = set()

    file_name = str(Path(root) / Path(file_str))

    # Check for correct file names and endings for each type, create file path, and return "key" to be used for column name in path df later with "value" being file path
    if "_rc" in file_name.lower() and Path(file_str).suffix == '.csv':
        key = 'count_file_path'

        rc_table_path = (Path(root) / file_str)

        rc_data = pd.read_csv(rc_table_path)
        print(f"Found RC file: {file_str}")

        # Make count table tab-delimited ending with ".count.txt" for mageck_test() recognition
        # TODO: move this to mageck_count() as it is currently making count tables for read count tables that have no cnttbl?
        count_file = file_str.replace(".csv", ".count.txt")
        rc_data.to_csv(
            (Path(root) / count_file),
            sep="\t",
            index=False,
        )
        value = str(Path(root) / count_file)
        print(f"Made count table for: {file_str}")

        return key, value
    
    elif "_cnttbl" in file_name.lower() and Path(file_str).suffix in ['.txt', '.tsv']:
        key = 'cnttbl_path'
        value = str(Path(root) / file_str)
        print(f"Found cnttbl file: {file_str}")
        return key, value

    elif "fastq" in file_name.lower() and Path(file_str).suffix in ['.fastq', '.fq']:
        key = 'fastq_dir_path'
        if root not in fastq_dirs:
            fastq_dirs.add(root)
            value = str(Path(root))
            print(f"Found fastq directory: {root}")
            return key, value

    elif Path(file_str).suffix == '.xlsx':
        print(f"Convert all .xlsx files to .csv before running")
        print(f"Skipping .xlsx file: {file_str}")
    elif '.count' in file_name.lower():
        pass
    else:
        raise FileNotFoundError("Make sure fastq file ends in .fastq or .fq, cnttbl is tab-delimited and ends in .txt or .txv, and read count table is comma-delimited ending with .csv not .xlsx.")
    
@profile(stdout=False, filename = r'C:\Users\spejo\Documents\2_CRISPR_analysis_test_output\crispr_analysis_pipeline_baseline.prof')
def run_pipeline(input_dir:str, output_dir:str, overwrite:bool = False):
    library = r"C:\Users\spejo\Documents\1_CRISPR_analysis_test_input\FASTQ\TKOv3_guide_sequences_key.csv"
    
    # Initialize DataFrames with columns
    rc_paths_df = pd.DataFrame(columns=['count_file_path', 'cnttbl_path', 'rc_out_dir'])
    fastq_paths_df = pd.DataFrame(columns=['fastq_dir_path', 'cnttbl_path', 'fastq_out_dir'])
    
    # Create list of screens that have already been analyzed
    analyzed_files = check_existing_files(input_dir, output_dir)
    
    # Track directories we've processed to avoid duplicates
    processed_dirs = set()
    
    for root, dirs, files in os.walk(input_dir):
        # Skip if no files in directory
        if not files:
            continue
            
        if "read_count" in root.lower():
            # Process read count directory
            if any('_cnttbl' in f.lower() for f in files) and any('_rc' in f.lower() and '.csv' in f.lower() for f in files):
                rc_files = []
                cnttbl = None
                rc_outpath = reverse_dir(input_dir=input_dir, root=root, output_dir=output_dir)
                
                for file in files:
                    if '_cnttbl' in file.lower():
                        cnttbl = str(Path(root) / file)
                    elif '_rc' in file.lower() and '.csv' in file.lower():
                        file_path = str(Path(root) / file)
                        # Skip if file has already been analyzed and overwrite is False
                        if file_path in analyzed_files and not overwrite:
                            logging.info(f"Skipping {file} - already analyzed")
                            continue
                        rc_files.append(file_path)
                
                if cnttbl and rc_files:
                    for rc_file in rc_files:
                        rc_paths_df.loc[len(rc_paths_df)] = {
                            'count_file_path': rc_file,
                            'cnttbl_path': cnttbl,
                            'rc_out_dir': str(rc_outpath)
                        }
                    logging.info(f"Added read count files from {root} for analysis")
                    
        elif "fastq" in root.lower():
            # Process FASTQ directory
            if root not in processed_dirs and any('_cnttbl' in f.lower() for f in files) and any(f.lower().endswith(('.fq', '.fastq')) for f in files):
                # Skip if directory has already been analyzed and overwrite is False
                if str(Path(root)) in analyzed_files and not overwrite:
                    logging.info(f"Skipping {root} - already analyzed")
                    continue
                
                cnttbl = None
                fastq_outpath = reverse_dir(input_dir=input_dir, root=root, output_dir=output_dir)
                
                # Find contrast table
                for file in files:
                    if '_cnttbl' in file.lower():
                        cnttbl = str(Path(root) / file)
                        break
                
                if cnttbl:
                    fastq_paths_df.loc[len(fastq_paths_df)] = {
                        'fastq_dir_path': root,
                        'cnttbl_path': cnttbl,
                        'fastq_out_dir': str(fastq_outpath)
                    }
                    processed_dirs.add(root)
                    logging.info(f"Added FASTQ directory {root} for analysis")

    # Process read count files
    if not rc_paths_df.empty:
        for _, row in rc_paths_df.iterrows():
            if not validate_input_tables(contrasts=row.cnttbl_path, counts=row.count_file_path):
                logging.warning(f"Skipping analysis for files in {row.cnttbl_path} due to validation error.")
                continue

            # Skip if file has already been analyzed and overwrite is False
            if row.count_file_path in analyzed_files and not overwrite:
                logging.info(f"Skipping analysis for {row.count_file_path} - already analyzed")
                continue

            test_files = mageck_test(input_file=row.count_file_path, 
                                   contrasts=row.cnttbl_path, 
                                   output_dir=row.rc_out_dir)
            logging.info(f"MAGeCK test results written to: {', '.join(test_files)}")
            
            drugz_files = run_drugz(input_file=row.count_file_path, 
                                  contrasts=row.cnttbl_path, 
                                  output_dir=row.rc_out_dir)
            logging.info(f"DrugZ results written to: {', '.join(drugz_files)}")
            logging.info(f"Finished analysis on read count files in {row.rc_out_dir}")

    # Process FASTQ files
    if not fastq_paths_df.empty:
        for _, row in fastq_paths_df.iterrows():
            # Skip if directory has already been analyzed and overwrite is False
            if row.fastq_dir_path in analyzed_files and not overwrite:
                logging.info(f"Skipping analysis for {row.fastq_dir_path} - already analyzed")
                continue

            sample_sheet = generate_sample_sheet(row.fastq_dir_path)
            logging.info(f"Generated sample sheet for {row.fastq_dir_path} with {len(sample_sheet)} samples")

            if not validate_input_tables(contrasts=row.cnttbl_path, sample_sheet=sample_sheet):
                logging.warning(f"Skipping analysis for files in {row.cnttbl_path} due to validation error.")
                continue

            count_file = mageck_count(
                sample_sheet=sample_sheet, 
                library_file=library,
                output_dir=row.fastq_out_dir,
            )
            logging.info(f"Count results written to: {count_file}")
            
            # Process the count file with all contrasts
            test_files = mageck_test(
                input_file=count_file,
                contrasts=row.cnttbl_path, 
                output_dir=row.fastq_out_dir
            )
            logging.info(f"Test results written to: {', '.join(test_files)}")
            
            drugz_files = run_drugz(
                input_file=count_file, 
                contrasts=row.cnttbl_path, 
                output_dir=row.fastq_out_dir
            )
            logging.info(f"DrugZ results written to: {', '.join(drugz_files)}")
            
            logging.info(f"Finished analysis on FASTQ files in {row.fastq_out_dir}")

    logging.info("Pipeline execution completed")


# TODO: could change mageck_count() and mageck_test() so count table output is easier to find/read

if __name__ == "__main__":
    # Example usage
    input_folder = r"C:\Users\spejo\Documents\1_CRISPR_analysis_test_input"
    output_folder = r"C:\Users\spejo\Documents\2_CRISPR_analysis_test_output"

    run_pipeline(input_folder, output_folder, overwrite=False)